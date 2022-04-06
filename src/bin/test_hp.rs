extern crate dbf_beam_simulator;
use scorus::healpix::utils::nside2npix;

use dbf_beam_simulator::utils::integrate_az;

fn main() {
    let nside = 16;
    let npix = nside2npix(nside);
    let map: Vec<_> = (0..npix).map(|_ipix| 1.0).collect();

    let (m, w, theta) = integrate_az(&map);

    for (&m1, (&w1, &theta1)) in m.iter().zip(w.iter().zip(theta.iter())) {
        println!("{} {} {}", m1, w1, theta1);
    }
}
