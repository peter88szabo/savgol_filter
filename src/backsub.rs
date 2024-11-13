fn lubksb(a: &[Vec<f64>], n: usize, _np: usize, indx: &[usize], b: &mut Vec<f64>) {
    let mut ii = 0;

    // Forward substitution
    for i in 0..n {
        let ll = indx[i];
        let mut sum = b[ll];
        b[ll] = b[i];
        
        if ii != 0 {
            for j in ii..i {
                sum -= a[i][j] * b[j];
            }
        } else if sum != 0.0 {
            ii = i + 1;
        }
        
        b[i] = sum;
    }

    // Back substitution
    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in (i + 1)..n {
            sum -= a[i][j] * b[j];
        }
        b[i] = sum / a[i][i];
    }
}

