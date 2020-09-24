package quasiNewton

import matrix.*
import java.lang.IllegalArgumentException
import kotlin.math.abs
import kotlin.math.sqrt

fun Matrix.BFSG_update(y: ColumnVector, s: ColumnVector): Unit {
    val sby = s - this * y
    val a = dot(y, s)
    this += (sby * s.transpose() + s * sby.transpose()) / a
    this -= (s * s.transpose()) * (dot(sby, y) / a / a)
}

class Pointwithfg(var x: ColumnVector, val f: (VecTor) -> Double, val grad: (VecTor) -> ColumnVector) {
    val F by lazy { f(x) }
    val Grad by lazy { grad(x) }

    fun SearchWolfe(delta: ColumnVector, c1: Double = 0.25, c2: Double = 0.75): Pair<Pointwithfg, ColumnVector> {
        var left = this
        var right = Pointwithfg(x + delta, f, grad)
        lateinit var middle: Pointwithfg
        do {
            if (right.F > F + dot(Grad, right.x - x) * c1 || (right.F >= left.F) && (left != this)) break
            if (abs(dot(right.Grad, delta)) <= -dot(Grad, delta) * c2) return Pair(right, right.x - x)
            if (dot(right.Grad, delta) >= 0) {
                val tmp = left;left = right;right = tmp;
                break
            }
            left = right
            right = Pointwithfg(right.x * 2.0 - x, f, grad)
        } while (true)
        while ((right.x - left.x).norm() > 0) {
            middle = Pointwithfg((left.x + right.x) / 2.0, f, grad)
            if (right.F > F + dot(Grad, middle.x - x) * c1 || middle.F >= left.F) {
                right = middle
            } else {
                if (abs(dot(middle.Grad, delta)) <= -dot(Grad, delta) * c2) return Pair(middle, middle.x - x)
                if (dot(middle.Grad, right.x - left.x) >= 0) right = left
                left = middle
            }
        }
        return Pair(middle, middle.x - x)
        /*var left: Double = 0.0
        var right: Double = 1.0
        var lamda: Double
        var fleft = F
        var xnew: ColumnVector
        var fnew: Double
        var gnew: ColumnVector
        do {
            xnew = x + delta * right
            fnew = f(xnew)
            if (fnew > F + dot(Grad, delta) * c1 * right || (fnew >= fleft) && (left > 0)) break
            gnew = grad(xnew)
            if (abs(dot(gnew, delta)) <= -dot(Grad, delta) * c2) {
                x = xnew
                F = fnew
                Grad = gnew
                return delta * right
            }
            if (dot(gnew, delta) >= 0) {
                fleft = fnew
                val tmp = left;left = right;right = tmp;
                break
            }
            left = right
            fleft = fnew
            right *= 2
        } while (true)
        while (left < right) {
            lamda = (left + right) / 2
            xnew = x + delta * lamda
            fnew = f(xnew)
            if (fnew > F + dot(Grad, delta) * c1 * lamda || fnew >= fleft) right = lamda
            else {
                gnew = grad(xnew)
                if (abs(dot(gnew, delta)) <= -dot(Grad, delta) * c2) {
                    x = xnew
                    F = fnew
                    Grad = gnew
                    return delta * lamda
                }
                if (dot(gnew, delta) * (right - left) >= 0) right = left
                left = lamda
                fleft = fnew
            }
        }*/
    }
}

fun fmin(f: (VecTor) -> Double, grad: (VecTor) -> ColumnVector, x0: ColumnVector): Pointwithfg {
    val d = x0.size
    val B = Identity(d)
    val err = sqrt(d.toDouble()) * 1e-6
    var point = Pointwithfg(x0, f, grad)

    // with(point) {
    while (point.Grad.norm() > err) {
        val delta = -B * point.Grad
        var OldGrad = point.Grad
        val (newpoint, s) = point.SearchWolfe(delta)
        point=newpoint
        B.BFSG_update(point.Grad - OldGrad, s)
    }
    //}
    return point
}