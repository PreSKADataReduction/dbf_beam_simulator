#![allow(unused_variables)]

extern crate dbf_beam_simulator as lpda;




fn main(){
    let nside=256;
    let (ant_x, ant_y):(Vec<_>, Vec<_>)=[(0.0, -2.0), (0.0, -1.0), (0.0, 0.0), (0.0, 1.0), (0.0, 2.0)].into_iter().unzip();
    let ant_z=vec![0.0; ant_x.len()];
    let ant_w=vec![1.0; ant_x.len()];
    //let array_beam=calc_array_beam(nside, &ant_x, &ant_y, &ant_z, &ant_w, 50e6);
    //write_map("array_beam.fits", &[&array_beam], false, true);
}