use num::{Complex, Zero};

use rustfft::{FftDirection, FftNum, FftPlanner};

use ndarray::{s, Array2, ArrayView2, ArrayViewMut, ArrayViewMut2, Dimension};

pub fn fftshift2<T>(in_data: ArrayView2<T>) -> Array2<T>
where
    T: Copy,
{
    assert!(in_data.shape()[0] % 2 == 0);
    let h = in_data.shape()[0];
    let w = in_data.shape()[1];
    let mut result =
        unsafe { Array2::uninit((in_data.shape()[0], in_data.shape()[1])).assume_init() };

    result
        .slice_mut(s![0..h / 2, 0..w / 2])
        .assign(&in_data.slice(s![h / 2..h, w / 2..w]));
    result
        .slice_mut(s![0..h / 2, w / 2..w])
        .assign(&in_data.slice(s![h / 2..h, 0..w / 2]));

    result
        .slice_mut(s![h / 2..h, 0..w / 2])
        .assign(&in_data.slice(s![0..h / 2, w / 2..w]));
    result
        .slice_mut(s![h / 2..h, w / 2..w])
        .assign(&in_data.slice(s![0..h / 2, 0..w / 2]));
    result
}

fn _fft<T: FftNum>(input: &mut [Complex<T>], output: &mut [Complex<T>], inverse: bool) {
    let mut planner = FftPlanner::new();
    let len = input.len();
    let fft = planner.plan_fft(
        len,
        if inverse {
            FftDirection::Inverse
        } else {
            FftDirection::Forward
        },
    );
    let scratch_len = fft.get_outofplace_scratch_len();
    let mut scratch = vec![Complex::zero(); scratch_len];
    fft.process_outofplace_with_scratch(input, output, &mut scratch);
}

pub fn fft<T: FftNum>(input: &mut [Complex<T>], output: &mut [Complex<T>]) {
    _fft(input, output, false);
}

pub fn ifft<T: FftNum>(input: &mut [Complex<T>], output: &mut [Complex<T>]) {
    _fft(input, output, true);
    /*
    for v in output.iter_mut() {
        *v = v.unscale(T::from(input.len() as u32));
    }*/
}

pub fn fft2(input: ArrayViewMut2<Complex<f64>>, output: ArrayViewMut2<Complex<f64>>) {
    fftnd(input, output, &[0, 1]);
}

pub fn ifft2(input: ArrayViewMut2<Complex<f64>>, output: ArrayViewMut2<Complex<f64>>) {
    ifftnd(input, output, &[1, 0]);
}

pub fn fftn<T: FftNum, D: Dimension>(
    input: &mut ArrayViewMut<Complex<T>, D>,
    output: &mut ArrayViewMut<Complex<T>, D>,
    axis: usize,
) {
    _fftn(input, output, axis, false);
}

pub fn ifftn<T: FftNum, D: Dimension>(
    input: &mut ArrayViewMut<Complex<T>, D>,
    output: &mut ArrayViewMut<Complex<T>, D>,
    axis: usize,
) {
    _fftn(input, output, axis, true);
}

fn _fftn<T: FftNum, D: Dimension>(
    input: &mut ArrayViewMut<Complex<T>, D>,
    output: &mut ArrayViewMut<Complex<T>, D>,
    axis: usize,
    inverse: bool,
) {
    if inverse {
        mutate_lane(input, output, ifft, axis)
    } else {
        mutate_lane(input, output, fft, axis)
    }
}

pub fn fftnd<T: FftNum, D: Dimension>(
    mut input: ArrayViewMut<Complex<T>, D>,
    mut output: ArrayViewMut<Complex<T>, D>,
    axes: &[usize],
) {
    _fftnd(&mut input, &mut output, axes, false);
}

pub fn ifftnd<T: FftNum, D: Dimension>(
    mut input: ArrayViewMut<Complex<T>, D>,
    mut output: ArrayViewMut<Complex<T>, D>,
    axes: &[usize],
) {
    _fftnd(&mut input, &mut output, axes, true);
}

fn _fftnd<T: FftNum, D: Dimension>(
    input: &mut ArrayViewMut<Complex<T>, D>,
    output: &mut ArrayViewMut<Complex<T>, D>,
    axes: &[usize],
    inverse: bool,
) {
    let len = axes.len();
    for (i, &axis) in axes.iter().enumerate() {
        _fftn(input, output, axis, inverse);
        if i < len - 1 {
            let mut outrows = output.rows_mut().into_iter();
            for mut row in input.rows_mut() {
                let mut outrow = outrows.next().unwrap();
                row.as_slice_mut()
                    .unwrap()
                    .copy_from_slice(outrow.as_slice_mut().unwrap());
            }
        }
    }
}

fn mutate_lane<T: Zero + Clone, D: Dimension>(
    input: &mut ArrayViewMut<T, D>,
    output: &mut ArrayViewMut<T, D>,
    f: fn(&mut [T], &mut [T]) -> (),
    axis: usize,
) {
    if axis > 0 {
        input.swap_axes(0, axis);
        output.swap_axes(0, axis);
        {
            let mut outrows = output.rows_mut().into_iter();
            for row in input.rows_mut() {
                let mut outrow = outrows.next().unwrap();
                let mut vec = row.to_vec();
                let mut out = vec![Zero::zero(); outrow.len()];
                f(&mut vec, &mut out);
                for i in 0..outrow.len() {
                    outrow[i] = out.remove(0);
                }
            }
        }
        input.swap_axes(0, axis);
        output.swap_axes(0, axis);
    } else {
        let mut outrows = output.rows_mut().into_iter();
        for mut row in input.rows_mut() {
            let mut outrow = outrows.next().unwrap();
            f(row.as_slice_mut().unwrap(), outrow.as_slice_mut().unwrap());
        }
    }
}
