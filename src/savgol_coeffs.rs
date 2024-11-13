// Maximum polynomial order
fn savgol(
    data: &[f64],    // Input data series to be smoothed
    nl: usize,       // Number of points leftward (past points)
    nr: usize,       // Number of points rightward (future points)
    ld: usize,       // Order of derivative (0 = smoothed function)
    m: usize         // Order of the smoothing polynomial (usual values: m = 2 or 4)
) -> Vec<f64> {
    let np = data.len(); // Length of the data series

    if np < nl + nr + 1 || nl < 0 || nr < 0 || ld > m | nl + nr < m {
        panic!("Invalid arguments in savgol");
    }

    // Initialize variables
    let mut a = vec![vec![0.0; m + 1]; MMAX + 1];
    let mut b = vec![0.0; m + 1];
    let mut indx = vec![0; m + 1];
    let mut c = vec![0.0; np]; // Initialize the smoothed data vector

    // Construct the smoothing polynomial matrix
    for ipj in 0..=2 * m {
        let mut sum = if ipj == 0 { 1.0 } else { 0.0 };

        for k in 1..=nr {
            sum += (k as f64).powi(ipj as i32);
        }
        for k in 1..=nl {
            sum += (-k as f64).powi(ipj as i32);
        }

        let mm = ipj.min(2 * m - ipj);
        for imj in (-mm..=mm).step_by(2) {
            a[(ipj + imj) / 2][(ipj - imj) / 2] = sum;
        }
    }

    let d = ludcmp(&mut a, m + 1, m + 1, &mut indx);

    for j in 0..=m {
        b[j] = 0.0;
    }
    b[ld] = 1.0;

    lubksb(&a, m + 1, m + 1, &indx, &mut b);

    // Apply smoothing coefficients to data series
    for kk in 0..np {
        c[kk] = 0.0;
        for k in -(nl as i32)..=nr as i32 {
            let mut sum = b[0];
            let mut fac = 1.0;

            for mm in 1..=m {
                fac *= k as f64;
                sum += b[mm] * fac;
            }

            let idx = (kk as i32 - k).rem_euclid(np as i32) as usize;
            c[kk] += sum * data[idx];
        }
    }

    c // the smoothed data vector
}

