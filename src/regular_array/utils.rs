use ndarray::{s, Array1, Array2, ArrayView2};

pub fn full2quarter<T>(full: ArrayView2<T>) -> Array2<T>
where
    T: Copy,
{
    let h = full.shape()[0];
    let w = full.shape()[1];
    full.slice(s![h / 2..h, w / 2..w]).to_owned()
}

pub fn quarter2full<T>(quarter: ArrayView2<T>) -> Array2<T>
where
    T: Copy + Default,
{
    let h = quarter.shape()[0] * 2;
    let w = quarter.shape()[1] * 2;
    let mut full = Array2::<T>::default((h, w));
    full.slice_mut(s![h / 2..h, w / 2..w]).assign(&quarter);

    
    full.slice_mut(s![1..h/2;-1, 1..w/2;-1])
        .assign(&quarter.slice(s![1..h / 2, 1..w / 2]));
    full.slice_mut(s![h/2..h, 1..w/2;-1])
        .assign(&quarter.slice(s![.., 1..w / 2]));
    full.slice_mut(s![1..h/2;-1, w/2..w])
        .assign(&quarter.slice(s![1..h / 2, ..]));
    full.slice(s![1.., 1..]).to_owned()

    
}

pub fn flattern_quarter_wgt(wgt: ArrayView2<f64>) -> Vec<f64> {
    wgt.iter().skip(1).cloned().collect()
}

pub fn deflattern_quarter_wgt(wgt: &[f64], h: usize, w: usize) -> Array2<f64> {
    let one = [1.0];
    if wgt.len()+1!=(h+1)/2*((w+1)/2)
    {
        println!("{} {} {}", wgt.len()+1, h, w);
    }
    Array1::from_iter(one.iter().chain(wgt.iter()).cloned())
        .into_shape(((h+1) / 2, (w+1) / 2))
        .unwrap()
}
