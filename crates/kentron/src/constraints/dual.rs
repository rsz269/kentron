//! Dual numbers: a value with its derivatives along `N` directions, carried through
//! arithmetic exactly (forward-mode differentiation), so a residual written once gives its
//! value and its Jacobian row together.

use std::ops::{Add, Div, Mul, Neg, Sub};

/// A value and its partial derivatives with respect to `N` variables.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Dual<const N: usize> {
    pub value: f64,
    pub derivatives: [f64; N],
}

impl<const N: usize> Dual<N> {
    /// A constant: no derivatives.
    pub fn constant(value: f64) -> Self {
        Dual { value, derivatives: [0.0; N] }
    }

    /// Variable `index` of the `N`, at `value`.
    pub fn variable(value: f64, index: usize) -> Self {
        let mut derivatives = [0.0; N];
        derivatives[index] = 1.0;
        Dual { value, derivatives }
    }

    /// `f(self)` given `f`'s value and derivative at `self.value`.
    fn chain(self, value: f64, slope: f64) -> Self {
        let mut derivatives = self.derivatives;
        for d in &mut derivatives {
            *d *= slope;
        }
        Dual { value, derivatives }
    }

    pub fn sin(self) -> Self {
        self.chain(self.value.sin(), self.value.cos())
    }

    pub fn cos(self) -> Self {
        self.chain(self.value.cos(), -self.value.sin())
    }

    /// The square root; its derivative at 0 taken as 0 (where a distance vanishes).
    pub fn sqrt(self) -> Self {
        let root = self.value.sqrt();
        self.chain(root, if root > 0.0 { 0.5 / root } else { 0.0 })
    }

    /// The absolute value; its derivative at 0 taken as 0.
    pub fn abs(self) -> Self {
        let sign = if self.value > 0.0 { 1.0 } else if self.value < 0.0 { -1.0 } else { 0.0 };
        self.chain(self.value.abs(), sign)
    }

    /// `atan2(self, x)`: the angle of the vector (x, self).
    pub fn atan2(self, x: Self) -> Self {
        let y = self;
        let r2 = x.value * x.value + y.value * y.value;
        let mut derivatives = [0.0; N];
        if r2 > 0.0 {
            for (i, d) in derivatives.iter_mut().enumerate() {
                *d = (x.value * y.derivatives[i] - y.value * x.derivatives[i]) / r2;
            }
        }
        Dual { value: y.value.atan2(x.value), derivatives }
    }
}

impl<const N: usize> Add for Dual<N> {
    type Output = Self;
    fn add(mut self, other: Self) -> Self {
        self.value += other.value;
        for (a, b) in self.derivatives.iter_mut().zip(other.derivatives) {
            *a += b;
        }
        self
    }
}

impl<const N: usize> Sub for Dual<N> {
    type Output = Self;
    fn sub(mut self, other: Self) -> Self {
        self.value -= other.value;
        for (a, b) in self.derivatives.iter_mut().zip(other.derivatives) {
            *a -= b;
        }
        self
    }
}

impl<const N: usize> Mul for Dual<N> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let mut derivatives = [0.0; N];
        for (i, d) in derivatives.iter_mut().enumerate() {
            *d = self.derivatives[i] * other.value + self.value * other.derivatives[i];
        }
        Dual { value: self.value * other.value, derivatives }
    }
}

impl<const N: usize> Div for Dual<N> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let inverse = 1.0 / other.value;
        let value = self.value * inverse;
        let mut derivatives = [0.0; N];
        for (i, d) in derivatives.iter_mut().enumerate() {
            *d = (self.derivatives[i] - value * other.derivatives[i]) * inverse;
        }
        Dual { value, derivatives }
    }
}

impl<const N: usize> Neg for Dual<N> {
    type Output = Self;
    fn neg(mut self) -> Self {
        self.value = -self.value;
        for d in &mut self.derivatives {
            *d = -*d;
        }
        self
    }
}

impl<const N: usize> Add<f64> for Dual<N> {
    type Output = Self;
    fn add(mut self, other: f64) -> Self {
        self.value += other;
        self
    }
}

impl<const N: usize> Sub<f64> for Dual<N> {
    type Output = Self;
    fn sub(mut self, other: f64) -> Self {
        self.value -= other;
        self
    }
}

impl<const N: usize> Mul<f64> for Dual<N> {
    type Output = Self;
    fn mul(mut self, other: f64) -> Self {
        self.value *= other;
        for d in &mut self.derivatives {
            *d *= other;
        }
        self
    }
}

impl<const N: usize> Div<f64> for Dual<N> {
    type Output = Self;
    fn div(self, other: f64) -> Self {
        self * (1.0 / other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The derivatives of x·y / (1 + x²) + sin(x) at (0.5, 2) by hand: with respect to x,
    /// y(1 - x²)/(1 + x²)² + cos x; with respect to y, x/(1 + x²).
    #[test]
    fn derivatives_of_a_rational_and_trigonometric_expression() {
        let (x, y) = (Dual::<2>::variable(0.5, 0), Dual::<2>::variable(2.0, 1));
        let f = x * y / (x * x + 1.0) + x.sin();
        let (xv, yv) = (0.5_f64, 2.0_f64);
        let dx = yv * (1.0 - xv * xv) / (1.0 + xv * xv).powi(2) + xv.cos();
        let dy = xv / (1.0 + xv * xv);
        assert!((f.derivatives[0] - dx).abs() < 1e-14 && (f.derivatives[1] - dy).abs() < 1e-14, "{f:?}");
    }

    /// atan2's derivatives: the angle of (x, y) turns at 1/r per unit across the radius.
    #[test]
    fn derivatives_of_an_angle() {
        let (y, x) = (Dual::<2>::variable(1.0, 0), Dual::<2>::variable(1.0, 1));
        let a = y.atan2(x);
        assert!((a.value - std::f64::consts::FRAC_PI_4).abs() < 1e-15);
        assert!((a.derivatives[0] - 0.5).abs() < 1e-15 && (a.derivatives[1] + 0.5).abs() < 1e-15, "{a:?}");
    }
}

/// A value with its gradient and Hessian along `N` variables (second-order forward-mode
/// differentiation), for Newton's method.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Jet<const N: usize> {
    pub value: f64,
    pub gradient: [f64; N],
    /// Symmetric, stored whole.
    pub hessian: [[f64; N]; N],
}

impl<const N: usize> Jet<N> {
    pub fn constant(value: f64) -> Self {
        Jet { value, gradient: [0.0; N], hessian: [[0.0; N]; N] }
    }

    pub fn variable(value: f64, index: usize) -> Self {
        let mut j = Jet::constant(value);
        j.gradient[index] = 1.0;
        j
    }

    /// 1 / self.
    fn reciprocal(self) -> Self {
        let inverse = 1.0 / self.value;
        let (i2, i3) = (inverse * inverse, inverse * inverse * inverse);
        let mut out = Jet::constant(inverse);
        for a in 0..N {
            out.gradient[a] = -self.gradient[a] * i2;
            for b in 0..N {
                out.hessian[a][b] = -self.hessian[a][b] * i2 + 2.0 * self.gradient[a] * self.gradient[b] * i3;
            }
        }
        out
    }
}

impl<const N: usize> Add for Jet<N> {
    type Output = Self;
    fn add(mut self, o: Self) -> Self {
        self.value += o.value;
        for a in 0..N {
            self.gradient[a] += o.gradient[a];
            for b in 0..N {
                self.hessian[a][b] += o.hessian[a][b];
            }
        }
        self
    }
}

impl<const N: usize> Sub for Jet<N> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self + (-o)
    }
}

impl<const N: usize> Neg for Jet<N> {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}

impl<const N: usize> Mul for Jet<N> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let mut out = Jet::constant(self.value * o.value);
        for a in 0..N {
            out.gradient[a] = self.gradient[a] * o.value + self.value * o.gradient[a];
            for b in 0..N {
                out.hessian[a][b] = self.hessian[a][b] * o.value + self.value * o.hessian[a][b] + self.gradient[a] * o.gradient[b] + self.gradient[b] * o.gradient[a];
            }
        }
        out
    }
}

impl<const N: usize> Div for Jet<N> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        self * o.reciprocal()
    }
}

impl<const N: usize> Add<f64> for Jet<N> {
    type Output = Self;
    fn add(mut self, o: f64) -> Self {
        self.value += o;
        self
    }
}

impl<const N: usize> Sub<f64> for Jet<N> {
    type Output = Self;
    fn sub(mut self, o: f64) -> Self {
        self.value -= o;
        self
    }
}

impl<const N: usize> Mul<f64> for Jet<N> {
    type Output = Self;
    fn mul(mut self, o: f64) -> Self {
        self.value *= o;
        for a in 0..N {
            self.gradient[a] *= o;
            for b in 0..N {
                self.hessian[a][b] *= o;
            }
        }
        self
    }
}

impl<const N: usize> Div<f64> for Jet<N> {
    type Output = Self;
    fn div(self, o: f64) -> Self {
        self * (1.0 / o)
    }
}

#[cfg(test)]
mod jet_tests {
    use super::*;

    /// f = x·y / (1 + x²) at (0.5, 2): its second derivatives by hand,
    /// f_xx = 2xy(x² − 3)/(1 + x²)³, f_xy = (1 − x²)/(1 + x²)², f_yy = 0.
    #[test]
    fn second_derivatives_of_a_rational_expression() {
        let (x, y) = (Jet::<2>::variable(0.5, 0), Jet::<2>::variable(2.0, 1));
        let f = x * y / (x * x + 1.0);
        let (xv, yv) = (0.5_f64, 2.0_f64);
        let fxx = 2.0 * xv * yv * (xv * xv - 3.0) / (1.0 + xv * xv).powi(3);
        let fxy = (1.0 - xv * xv) / (1.0 + xv * xv).powi(2);
        assert!((f.hessian[0][0] - fxx).abs() < 1e-14, "{f:?}");
        assert!((f.hessian[0][1] - fxy).abs() < 1e-14 && (f.hessian[1][0] - fxy).abs() < 1e-14 && f.hessian[1][1].abs() < 1e-15, "{f:?}");
    }
}
