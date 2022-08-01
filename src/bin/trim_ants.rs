use clap::{Arg, Command};

use std::{collections::BTreeMap, fs::File};

use serde_yaml::to_writer;

use serde::Serialize;

#[derive(Serialize)]
struct GridArrayCfg{
    pub full_array: Vec<(isize, isize)>,
    pub full_bl: Vec<(isize, isize, usize)>,
    pub trimed_array: Vec<(isize, isize)>,
    pub trimed_bl: Vec<(isize, isize, usize)>
}

pub fn regulate_baseline(bl: (isize, isize)) -> (isize, isize) {
    /*
    if bl.0 > 0 || (bl.0==0 && bl.1 >0){
        bl
    } else {
        (-bl.0, -bl.1)
    }*/
    bl
}

pub fn get_baseline(ants: &[(isize, isize)]) -> BTreeMap<(isize, isize), usize> {
    let mut result = BTreeMap::new();
    for a1 in ants {
        for a2 in ants {
            let bl0 = (a2.0 - a1.0, a2.1 - a1.1);
            let bl = regulate_baseline(bl0);

            *result.entry(bl).or_insert(0) += 1;
        }
    }
    
    result
}

pub fn get_ants(array_size: usize) -> Vec<(isize, isize)> {
    let array_size = array_size as isize;
    //assert!(array_size % 2 == 0);
    let mut result = Vec::new();
    for i in (-(array_size + 0) / 2)..(array_size + 1) / 2 {
        //println!("{}", i);
        for j in (-(array_size + 0) / 2)..(array_size + 1) / 2 {
            result.push((i, j));
        }
    }
    result
}

pub fn sort_ants_by_grade(ants: &[(isize, isize)]) -> Vec<(isize, isize)> {
    let baseline = get_baseline(ants);
    let mut result = vec![0; ants.len()];
    for (i, a1) in ants.iter().enumerate() {
        for (j, a2) in ants.iter().enumerate() {
            let bl = regulate_baseline((a2.0 - a1.0, a2.1 - a1.1));
            if bl.0 != 0 || bl.1 != 0 {
                let w = baseline[&bl];
                result[i] += w;
                result[j] += w;
            }
        }
    }
    let mut result = ants
        .iter()
        .zip(result.iter())
        .map(|(&a, &w)| (a, w / 2))
        .collect::<Vec<_>>();
    result.sort_by_key(|x| (x.1 as isize));
    result.iter().cloned().map(|(a, _w)| a).collect()
}

pub fn trim_ants(ants: &[(isize, isize)]) -> Vec<(isize, isize)> {
    //let nants = ants.len();
    let mut sorted_ants = sort_ants_by_grade(ants);
    let bl = get_baseline(ants);
    loop {
        let bl1 = get_baseline(&sorted_ants[..sorted_ants.len() - 1]);
        if bl1.len() != bl.len() {
            break sorted_ants;
        }
        sorted_ants.pop();
    }
}

fn main() {
    let matches = Command::new("trim_ants")
        .arg(
            Arg::new("array_size")
                .short('s')
                .long("size")
                .takes_value(true)
                .required(true)
                .value_name("array size"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("out")
                .takes_value(true)
                .required(true)
                .value_name("outfile"),
        )
        .get_matches();

    //let mut cfg = BTreeMap::<&str, Vec<(isize, isize)>>::new();

    let mut outfile = File::create(matches.value_of("outfile").unwrap()).unwrap();
    let array_size = matches
        .value_of("array_size")
        .unwrap()
        .parse::<usize>()
        .unwrap();

    let ants = get_ants(array_size);
    //println!("{:?}", ants);
    

    let bl = get_baseline(&ants);
    let bl = bl.iter().map(|(&k, &v)| (k.0, k.1, v)).collect();
    

    let trimed = trim_ants(&ants);
    

    let trimed_bl = get_baseline(&trimed);
    let trimed_bl = trimed_bl.iter().map(|(&k, &v)| (k.0, k.1, v)).collect();
    

    let cfg=GridArrayCfg{
        full_array: ants,
        full_bl: bl,
        trimed_array: trimed,
        trimed_bl: trimed_bl
    };

    to_writer(&mut outfile, &cfg).unwrap();
}
