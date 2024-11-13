fn savgol_smooth(
    ndata: usize,
    signal: &mut [f64],
    mpoly: usize,
    nl: usize,
    nr: usize,
    nder: usize,
    hh: f64,
) {
    let nmax = 50;
    let mut coeff = vec![0.0; nmax];
    let mut ind = vec![0; nmax];
    let ysave = signal.to_vec(); // Save the original signal

    // Shift indices for given case `nl`, `nr`, `m`
    ind[0] = 0;

    // Indices for left shifts
    let mut j = 3;
    for i in 1..=nl {
        ind[i] = i as isize - j;
        j += 2;
    }

    // Indices for right shifts
    j = 2;
    for i in (nl + 1)..=(nl + nr) {
        ind[i] = (i as isize) - j;
        j += 2;
    }

    // Get the coefficients with Savitzky-Golay filter
    let mint = nl + nr + 1;
    savgol(&mut coeff, mint, nl, nr, nder, mpoly);

    // Apply filter to input data
    for i in 0..(ndata - nr) {
        signal[i] = 0.0;
        for j in 0..(nl + nr + 1) {
            let index = i as isize + ind[j];
            if index >= 0 {
                signal[i] += coeff[j] * ysave[index as usize];
            }
        }
    }

    // Adjust for derivative if needed
    if nder == 1 {
        for i in 0..(ndata - nr) {
            signal[i] /= hh;
        }
    } else if nder == 2 {
        for i in 0..(ndata - nr) {
            signal[i] *= 2.0 / (hh * hh);
        }
    }
}


