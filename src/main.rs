#![allow(unused_imports)]

use serde_yaml::to_writer;

use std::{
    fs::File
    , io::Write
};

use dbf_beam_simulator::{
    array_cfg::{
        AntCfg
        , ArrayCfg
    }
};


pub fn main() {
    let ant_cfg=AntCfg{
        pos:(0.0, 0.0, 0.0)
    };
    let array_cfg=ArrayCfg{
        ants:vec![ant_cfg.clone(), ant_cfg.clone()]
    };

    let mut outfile=File::create("a.yaml").unwrap();
    to_writer(&mut outfile, &array_cfg).unwrap();
}
