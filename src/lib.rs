#![allow(non_snake_case)]
use std::{
    f64::{
        consts::{
            PI
        }
    }
};

use serde::{Serialize, Deserialize};


use num::{
    complex::{
        Complex
    }
};

const LIGHT_SPEED:f64=2.99792458E8;

use scorus::{
    healpix::{
        utils::{
            nside2npix
            , npix2nside
            , nside2nring
            , nring2nside
        }
        , pix::{
            pix2vec_ring, pix2ring_ring
            , ring2z_ring
        }
        , rotation::rotate_ring
    }
    , coordinates::{
        SphCoord
        ,Vec3d
        , rotation3d::{
            RotMatrix
        }
    }
};

#[derive(Clone, Serialize, Deserialize)]
pub struct ArrayCfg{
    pub ants: Vec<AntCfg>
}

#[derive(Clone, Serialize, Deserialize)]
pub struct AntCfg{
    pub pos: (f64,f64,f64)
    , pub weight: f64
}

pub fn calc_array_beam(nside: usize, x_list: &[f64], y_list: &[f64], z_list: &[f64], w_list: &[f64], phi_list:&[f64], freq_Hz: f64, ground_cut: bool)->Vec<f64>{
    let npix=nside2npix(nside);
    let lambda=LIGHT_SPEED/(freq_Hz);
    (0..npix).map(|i|{
        if i<npix/2 || !ground_cut{
            let pointing=pix2vec_ring::<f64>(nside, i);
            x_list.iter().zip(y_list.iter().zip(z_list.iter().zip(w_list.iter().zip(phi_list.iter())))).map(|(&x, (&y, (&z, (&w, &phi)))) |{
                let dl=pointing[0]*x+pointing[1]*y+pointing[2]*z;
                let phase=dl/lambda*2.0*PI;

                Complex::from_polar(w, phase-phi)
            }).sum::<Complex<f64>>().norm_sqr()
        }else{
            0.0
        }
    }).collect()
}

pub fn calc_phase_from_pointing(x_list: &[f64], y_list: &[f64], z_list: &[f64], az_from_north: f64, zenith: f64, freq_Hz: f64)->Vec<f64>{
    let az_from_x=PI/2.0-az_from_north;
    let lambda=LIGHT_SPEED/freq_Hz;
    let dir=Vec3d::from_sph_coord(SphCoord::new(zenith, az_from_x));
    x_list.iter().zip(y_list.iter().zip(z_list.iter())).map(|(&x, (&y, &z))|{
        let dl=dir[0]*x+dir[1]*y+dir[2]*z;
        dl/lambda*2.0*PI
    }).collect()
}

pub fn integrate_az(hmap: &[f64])->(Vec<f64>, Vec<usize>, Vec<f64>){
    let npix=hmap.len();
    let nside=npix2nside(npix);
    let nring=nside2nring(nside);
    let mut wgt=vec![0_usize; nring];
    let mut mean_values=vec![0.0; nring];
    let ring_idx:Vec<_>=(0..npix).map(|ipix| pix2ring_ring(nside, ipix)).collect();
    ring_idx.iter().for_each(|&i| wgt[i-1]+=1);
    hmap.iter().zip(ring_idx.iter()).for_each(|(&x, &iring)|{
        mean_values[iring-1]+=x;
    });
    mean_values.iter_mut().zip(wgt.iter()).for_each(|(x, &w)|{
        *x/=w as f64
    });

    let theta=(1..=nring).map(|iring| ring2z_ring::<f64>(nside, iring).acos()).collect();
    (mean_values, wgt, theta)
}

pub fn calc_averaged_array_beam(lat_deg: f64, ant_beam: &[f64], x_list: &[f64], y_list: &[f64], z_list: &[f64], w_list: &[f64], phi_list:&[f64], freq_hz: f64)->(Vec<f64>, Vec<usize>, Vec<f64>)
{
    let nside=npix2nside(ant_beam.len());
    let array_beam=calc_array_beam(nside, x_list, y_list, z_list, w_list, phi_list, freq_hz, true);
    let total_beam:Vec<f64>=array_beam.iter().zip(ant_beam.iter()).map(|(&a, &b)|{a*b}).collect();
    let rot=RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 1.0, 0.0), (90.0-lat_deg).to_radians())
    *RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 0.0, 1.0), -90_f64.to_radians());
    let rotated_beam=rotate_ring(&total_beam, &rot);
    let (mean_beam, weight, theta)=integrate_az(&rotated_beam);
    (mean_beam, weight, theta)
}

pub fn averaged_beam_to_healpix(beam: &[f64])->Vec<f64>{
    let nring=beam.len();
    let nside=nring2nside(nring);
    let npix=nside2npix(nside);
    (0..npix).map(|ipix|{
        beam[pix2ring_ring(nside, ipix)-1]
    }).collect()
}

pub fn beam_opt_func_obj1(beam1: &[f64], beam2: &[f64], wgt: &[f64])->f64{
    assert_eq!(beam1.len(), wgt.len());
    assert_eq!(beam2.len(), wgt.len());
    let s1=beam1.iter().zip(wgt.iter()).map(|(&a,&b)|{a*b}).sum::<f64>();
    let s2=beam2.iter().zip(wgt.iter()).map(|(&a,&b)|{a*b}).sum::<f64>();
    beam1.iter().zip(beam2.iter().zip(wgt.iter())).map(|(&b1, (&b2, &w))|{
        (b1/s1-b2/s2)*w
    }).map(|x|x.powi(2)).sum::<f64>()
}

pub fn beam_opt_func_obj2(beam0: &[f64], 
    lat_deg: f64, ant_beam: &[f64], 
    x_list: &[f64], y_list: &[f64], z_list: &[f64], w_list: &[f64], phi_list:&[f64], freq_hz: f64)->f64
{
    let (mean_beam, weight, _theta)=calc_averaged_array_beam(lat_deg, ant_beam, x_list, y_list, z_list, w_list, phi_list, freq_hz);
    let weight:Vec<_>=weight.into_iter().map(|x| x as f64).collect();
    beam_opt_func_obj1(&beam0, &mean_beam, &weight)
}

pub fn zenith_ns_sym_array(y: &[f64], w: &[f64])->(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>){
    let mut ant_x=vec![0.0];
    let mut ant_y=vec![0.0];
    let mut ant_z=vec![0.0];
    let mut ant_w=vec![1.0];

    for (&y1, &w1) in y.iter().zip(w.iter()){
        ant_x.push(0.0);
        ant_y.push(y1);
        ant_z.push(0.0);
        ant_w.push(w1);

        ant_x.push(0.0);
        ant_y.push(-y1);
        ant_z.push(0.0);
        ant_w.push(w1);
    }
    (ant_x, ant_y, ant_z, ant_w)
}

pub fn sym_weight(w: &[f64])->Vec<f64>{
    let mut ant_w=vec![1.0];

    for &w1 in w{
        ant_w.push(w1);
        ant_w.push(w1);
    }
    ant_w
}

