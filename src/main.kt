import matrix.*
import matrix.Vector
import quasiNewton.fmin
import java.util.*

val n = 10
val A = ColumnVector(DoubleArray(n) { Random().nextDouble() })

fun simple(x: Vector) = dot(x, x) + dot(A, x)
fun hhh(x: Vector) = ColumnVector(DoubleArray(n){ i->x[i]*2+A[i]})
fun main() {
    val ans=fmin(::simple,::hhh,ColumnVector(DoubleArray(n) { Random().nextDouble() }))
    println(ans.F)
    println(ans.x.transpose())
}
