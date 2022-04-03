#![allow(non_snake_case)]

pub mod fft;

use std::{
    f64::{
        consts::{
            PI
        }
    }
};


use num::{
    complex::{
        Complex
    }
};

use ndarray::{
    Array1
    , Array2
    , ArrayView2
    , s
};

use crate::{
    fft::{
        fftshift2
        , fft2
    }
};

pub const LIGHT_SPEED:f64=2.99792458E8;

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
        , interp::{
            get_interpol_ring
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
    let npix=beam1.len();
    let domega=4.0*PI/npix as f64;
    let s1=beam1.iter().zip(wgt.iter()).map(|(&a,&b)|{a*b}).sum::<f64>();
    let s2=beam2.iter().zip(wgt.iter()).map(|(&a,&b)|{a*b}).sum::<f64>();
    beam1.iter().zip(beam2.iter().zip(wgt.iter())).map(|(&b1, (&b2, &w))|{
        (b1/s1-b2/s2)*w
    }).map(|x|x.powi(2)).sum::<f64>()/domega
}

pub fn beam_opt_func_obj2(beam0: &[f64], 
    lat_deg: f64, ant_beam: &[f64], 
    x_list: &[f64], y_list: &[f64], z_list: &[f64], w_list: &[f64], phi_list:&[f64], freq_hz: f64)->f64
{
    let (mean_beam, weight, _theta)=calc_averaged_array_beam(lat_deg, ant_beam, x_list, y_list, z_list, w_list, phi_list, freq_hz);
    let weight:Vec<_>=weight.into_iter().map(|x| x as f64).collect();
    beam_opt_func_obj1(&beam0, &mean_beam, &weight)
}

pub fn calc_averaged_ant_output(beam: &[f64], sky: &[f64])->f64{
    let beam_hp=averaged_beam_to_healpix(beam);
    assert_eq!(beam_hp.len(), sky.len());
    let norm=beam_hp.iter().sum::<f64>();
    beam_hp.iter().zip(sky.iter()).map(|(&b,&s)|{
        b*s
    }).sum::<f64>()/norm
}

pub fn calc_averaged_ant_output2(sky: &[f64], lat_deg: f64, ant_beam: &[f64],
    x_list: &[f64], y_list: &[f64], z_list: &[f64], w_list: &[f64], phi_list:&[f64], freq_hz: f64)->f64
{
    let (mean_beam, _weight, _theta)=calc_averaged_array_beam(lat_deg, ant_beam, x_list, y_list, z_list, w_list, phi_list, freq_hz);
    calc_averaged_ant_output(&mean_beam, &sky)
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

pub fn design_square_array(spacing: f64, freq_mega_hz: f64, sigma_deg: f64, n: isize)->Array2<f64>{
    let center=n/2;
    let lbd=LIGHT_SPEED/(freq_mega_hz*1e6);
    let sigma=spacing/lbd*sigma_deg.to_radians().sin();
    let mut beam_pattern=Array2::<Complex<f64>>::zeros((n as usize,n as usize));
    for i in 0..n{
        let x=i-center;
        let fx=x as f64/n as f64;
        let bx=f64::exp(-fx.powi(2)/(4.0*sigma.powi(2))); //4 for sqrt(B)
        for j in 0..n{
            let y=j-center;
            let fy=y as f64/n as f64;
            let by=f64::exp(-fy.powi(2)/(4.0*sigma.powi(2))); //4 for sqrt(B)
            beam_pattern[(i as usize,j as usize)]=Complex::<f64>::from(bx*by);
        }
    }    
    let mut beam_pattern=fftshift2(beam_pattern.view());
    let mut wgt=Array2::<Complex<f64>>::zeros((n as usize,n as usize));
    fft2(beam_pattern.view_mut(), wgt.view_mut());
    let norm=wgt[(0,0)].re;
    let mut wgt=wgt.map(|x| x.re);
    wgt.iter_mut().for_each(|x| *x=*x/norm);
    println!("{}", wgt[(1,1)]);
    println!("{}", wgt[(0,0)]);
    fftshift2(wgt.view())
}

pub fn pattern2wgt(hp: &[f64], d: f64, freq_mhz: f64, n: isize)->Array2<f64>{
    let npix=hp.len();
    let nside=npix2nside(npix);
    let center=n/2;
    let lbd=LIGHT_SPEED/(freq_mhz*1e6);
    let u=d/lbd;
    let mut projected=Array2::<Complex<f64>>::zeros((n as usize,n as usize));
    
    for i in 0..n{
        let fx=(i as f64-center as f64)/n as f64;
        let nx=fx/u;
        for j in 0..n{
            let fy=(j as f64-center as f64)/n as f64;
            let ny=fy/u;
            let r2=nx.powi(2)+ny.powi(2);
            if r2>1.0{
                continue;
            }
            let nz=(1.0-r2).sqrt();

            let (p,w)=get_interpol_ring(nside, SphCoord::from_xyz(nx, ny, nz));
            projected[(i as usize, j as usize)]=p.iter().zip(w.iter()).map(|(&i1,&w1)|{
                w1*hp[i1]
            }).sum::<f64>().sqrt().into();
        }
    }
    let mut projected=fftshift2(projected.view());
    let mut wgt=Array2::<Complex<f64>>::zeros((n as usize,n as usize));    
    fft2(projected.view_mut(), wgt.view_mut());
    let norm=wgt[(0,0)].re;
    let mut wgt=wgt.map(|x| x.re);
    wgt.iter_mut().for_each(|x| *x=*x/norm);
    fftshift2(wgt.view())        
}

pub fn wgt2pattern(wgt: ArrayView2<f64>, d: f64, freq_mhz: f64, nside: usize)->Vec<f64>{
    let npix=nside2npix(nside);
    let lbd=LIGHT_SPEED/(freq_mhz*1e6);
    let u=d/lbd;
    (0..npix).map(|ipix| {
        if ipix<npix/2{
            let Vec3d{x:nx,y:ny, z:_}=pix2vec_ring::<f64>(nside, ipix);
            let mut p=Complex::<f64>::new(0.0, 0.0);
            for i in 0..wgt.shape()[0]{
                let m=(i as isize-wgt.shape()[0] as isize/2) as f64;
                for j in 0..wgt.shape()[1]{
                    let w=wgt[(i,j)];
                    let n=(j as isize-wgt.shape()[1] as isize/2) as f64;
                    p+=Complex::<f64>::from_polar(w, 2.0*PI*(m*u*nx+n*u*ny));
                }
            }
            p.norm_sqr()
        }else{
            0.0
        }
    }).collect()
}

pub fn quarter_wgt2pattern(quarter_wgt: ArrayView2<f64>, d: f64, freq_mhz: f64, nside: usize)->Vec<f64>{
    let npix=nside2npix(nside);
    let lbd=LIGHT_SPEED/(freq_mhz*1e6);
    let u=d/lbd;
    (0..npix).map(|ipix| {
        if ipix<npix/2{
            let Vec3d{x:nx,y:ny, z:_}=pix2vec_ring::<f64>(nside, ipix);
            let mut p=0.0;
            for i in 0..quarter_wgt.shape()[0]{
                for j in 0..quarter_wgt.shape()[1]{
                    let w=quarter_wgt[(i,j)];
                    p+=w*(2.0*PI*i as f64*u*nx).cos()*(2.0*PI*j as f64*u*ny).cos()*
                    if i==0{
                        1.0
                    }else{
                        2.0
                    }*
                    if j==0{
                        1.0
                    }else{
                        2.0
                    };
                }
            }
            p.powi(2)
        }else{
            0.0
        }
    }).collect()
}

pub fn full2quarter<T>(full: ArrayView2<T>)->Array2<T>
where T:Copy
{
    let h=full.shape()[0];
    let w=full.shape()[1];
    full.slice(s![h/2..h, w/2..w]).to_owned()
}

pub fn quarter2full<T>(quarter: ArrayView2<T>)->Array2<T>
where T:Copy+Default{
    let h=quarter.shape()[0]*2;
    let w=quarter.shape()[1]*2;
    let mut full=Array2::<T>::default((h,w));
    full.slice_mut(s![h/2..h, w/2..w]).assign(&quarter);
    full.slice_mut(s![1..h/2;-1, 1..w/2;-1]).assign(&quarter.slice(s![1..h/2, 1..w/2]));
    full.slice_mut(s![h/2..h, 1..w/2;-1]).assign(&quarter.slice(s![.., 1..w/2]));
    full.slice_mut(s![1..h/2;-1, w/2..w]).assign(&quarter.slice(s![1..h/2, ..]));
    full
}

pub fn flattern_quarter_wgt(wgt: ArrayView2<f64>)->Vec<f64>{
    wgt.iter().skip(1).cloned().collect()
}

pub fn deflattern_quarter_wgt(wgt: &[f64], h: usize, w: usize)->Array2<f64>{
    let one=[1.0];
    Array1::from_iter(one.iter().chain(wgt.iter()).cloned()).into_shape((h/2, w/2)).unwrap()
}
