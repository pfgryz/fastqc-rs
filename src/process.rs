use chrono::{DateTime, Local};
use itertools::Itertools;
use needletail::{parse_fastx_file, Sequence};
use rustc_hash::FxHashMap as HashMap;
use serde_json::json;
use serde_json::Value;
use std::error::Error;
use std::ffi::OsStr;
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::Path;
use tera::{self, Context, Tera};

const BASES: [char; 5] = ['A', 'C', 'G', 'T', 'N'];
#[allow(unused)]
const A: usize = 0;
const C: usize = 1;
const G: usize = 2;
const T: usize = 3;
const N: usize = 4;

fn quartiles(hist: &[usize]) -> [f32; 5] {
    let sum = hist.iter().sum::<usize>();
    assert!(sum != 0);
    if sum == 1 {
        let value = hist.iter().enumerate().fold(
            0_usize,
            |acc, (value, &count)| if count > 0 { value } else { acc },
        ) as f32;
        return [value, value, value, value, value];
    }
    let mut ret = [0_f32; 5];
    // compute the quartiles
    for (i, quantile) in [0.25, 0.5, 0.75].iter().enumerate() {
        let rank = quantile * (sum - 1) as f64;
        let rank_ = rank.floor();
        let delta = rank - rank_;
        let n = rank_ as usize + 1;
        let mut acc = 0;
        let mut lo = None;
        for (hi, &count) in hist.iter().enumerate().filter(|(_, &count)| count > 0) {
            if acc == n && lo.is_some() {
                let lo = lo.unwrap() as f64;
                ret[i + 1] = (lo + (hi as f64 - lo) * delta) as f32;
                break;
            } else if acc + count > n {
                ret[i + 1] = hi as f32;
                break;
            }
            acc += count;
            lo = Some(hi);
        }
    }
    // compute lower, upper fences
    // TODO(lhepler): the UI reports these as min/max, which is incorrect.
    // If that's what we want, we can return those values.
    let iqr = ret[3] - ret[1];
    ret[0] = ret[1] - 1.5 * iqr;
    ret[4] = ret[3] + 1.5 * iqr;
    ret
}

pub(crate) fn process<P: AsRef<Path> + AsRef<OsStr>>(
    filename: P,
    k: u8,
    summary: Option<P>,
) -> Result<(), Box<dyn Error>> {
    let mut base_quality_count = HashMap::default();

    let mut reader = parse_fastx_file(&filename).expect("Invalid path/file");
    let mut broken_read = false;

    // Gather data from every record
    while let Some(record) = reader.next() {
        if let Ok(seqrec) = record {
            if let Some(qualities) = seqrec.qual() {
                for (pos, &q) in qualities.iter().enumerate() {
                    let rec = base_quality_count
                        .entry(pos)
                        .or_insert_with(|| vec![0_usize; 94]);
                    rec[q as usize - 33] += 1;
                }
            }
        } else {
            broken_read = true;
        }
    }
   
    // Data for base quality per position
    let mut base_quality_warn = "pass";
    let mut base_per_pos_data = Vec::new();
    for (position, qualities) in base_quality_count {
        let (sum, len) = qualities
            .iter()
            .enumerate()
            .fold((0_usize, 0_usize), |(s, l), (q, c)| (s + q * c, l + c));
        let avg = sum as f64 / len as f64;
        let values = quartiles(&qualities);
        if values.get(2).unwrap() <= &20_f32 {
            base_quality_warn = "fail"
        } else if values.get(2).unwrap() <= &25_f32 && base_quality_warn != "fail" {
            base_quality_warn = "warn"
        }
        base_per_pos_data.push(json!({
        "pos": position,
        "average": avg,
        "upper": values.get(4).unwrap(),
        "lower": values.get(0).unwrap(),
        "q1": values.get(1).unwrap(),
        "q3": values.get(3).unwrap(),
        "median":values.get(2).unwrap(),
        }));
    }

    let mut qpp_specs: Value =
        serde_json::from_str(include_str!("report/quality_per_pos_specs.json"))?;
    qpp_specs["data"]["values"] = json!(base_per_pos_data);

    let plots = json!({
        "base sequence quality": {"short": "base", "specs": qpp_specs.to_string()},
    });

    let file = Path::new(&filename).file_name().unwrap().to_str().unwrap();
    let meta = json!({
        "file name": {"name": "file name", "value": file},
        "canonical": {"name": "canonical", "value": "True"},
    });

    let mut templates = Tera::default();
    templates.register_filter("embed_source", embed_source);
    templates.add_raw_template("report.html.tera", include_str!("report/report.html.tera"))?;
    let mut context = Context::new();
    context.insert("plots", &plots);
    context.insert("meta", &meta);
    let local: DateTime<Local> = Local::now();
    context.insert("time", &local.format("%a %b %e %T %Y").to_string());
    context.insert("version", &env!("CARGO_PKG_VERSION"));
    context.insert("invalid_reads", &broken_read);
    let html = templates.render("report.html.tera", &context)?;
    io::stdout().write_all(html.as_bytes())?;

    if let Some(path) = summary {
        let output_path = Path::new(&path);
        templates.add_raw_template(
            "fastqc_summary.txt.tera",
            include_str!("report/fastqc_summary.txt.tera"),
        )?;
        context.insert("filename", &file);
        context.insert("base_quality_warn", &base_quality_warn);
        let txt = templates.render("fastqc_summary.txt.tera", &context)?;
        let mut file = File::create(output_path.join("fastqc_data.txt"))?;
        file.write_all(txt.as_bytes())?;
    }
    Ok(())
}

fn embed_source(
    value: &tera::Value,
    _: &std::collections::HashMap<String, tera::Value>,
) -> tera::Result<tera::Value> {
    let url = tera::try_get_value!("upper", "value", String, value);
    let source = reqwest::get(&url).unwrap().text().unwrap();
    Ok(tera::to_value(source).unwrap())
}

#[cfg(test)]
mod test {
    use super::quartiles;
    #[test]
    fn test_quartiles1() {
        let v1 = [-49.5, 24.75, 49.5, 74.25, 148.5];
        let v2 = quartiles(&vec![1; 100]);
        assert!(v1 == v2);
    }
    #[test]
    fn test_quartiles2() {
        let v1 = [6.25, 25.0, 25.0, 37.5, 56.25];
        let mut h = [0_usize; 76];
        h[25] = 75;
        h[75] = 25;
        let v2 = quartiles(&h);
        assert!(v1 == v2);
    }
    #[test]
    fn test_quartiles3() {
        let v1 = [43.75, 62.5, 75.0, 75.0, 93.75];
        let mut h = [0_usize; 76];
        h[25] = 25;
        h[75] = 75;
        let v2 = quartiles(&h);
        assert!(v1 == v2);
    }
}
