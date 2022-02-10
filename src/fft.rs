use fftw::{
    plan::{
        C2CPlan64
        , C2CPlan
    }
    , types::{
        Flag, Sign
        , c64
    }
};


use ndarray::{
    ArrayViewMut2
    , ArrayView2
    , Array2
    , s
};


pub fn fft2(mut in_data: ArrayViewMut2<c64>, 
    mut out_data: ArrayViewMut2<c64>){
    let mut plan=C2CPlan64::new(&[in_data.nrows(),in_data.ncols()], in_data.as_slice_mut().unwrap(), out_data.as_slice_mut().unwrap(), 
    Sign::Forward, Flag::ESTIMATE).unwrap();
    plan.c2c(&mut in_data.as_slice_mut().unwrap(), &mut out_data.as_slice_mut().unwrap()).unwrap();
}

pub fn ifft2(mut in_data: ArrayViewMut2<c64>, 
    mut out_data: ArrayViewMut2<c64>){
    let mut plan=C2CPlan64::new(&[in_data.nrows(),in_data.ncols()], in_data.as_slice_mut().unwrap(), out_data.as_slice_mut().unwrap(), 
    Sign::Backward, Flag::ESTIMATE).unwrap();
    plan.c2c(&mut in_data.as_slice_mut().unwrap(), &mut out_data.as_slice_mut().unwrap()).unwrap();    
}

pub fn fftshift2<T>(in_data: ArrayView2<T>) -> Array2<T>
where
    T: Copy,
{
    assert!(in_data.shape()[0] % 2 == 0);
    let h=in_data.shape()[0];
    let w=in_data.shape()[1];
    let mut result =
        unsafe { Array2::uninit((in_data.shape()[0], in_data.shape()[1])).assume_init() };

    result.slice_mut(s![0..h/2,0..w/2]).assign(&in_data.slice(s![h/2..h,w/2..w]));
    result.slice_mut(s![0..h/2, w/2..w]).assign(&in_data.slice(s![h/2..h, 0..w/2]));

    result.slice_mut(s![h/2..h,0..w/2]).assign(&in_data.slice(s![0..h/2,w/2..w]));
    result.slice_mut(s![h/2..h, w/2..w]).assign(&in_data.slice(s![0..h/2, 0..w/2]));
    result
}
