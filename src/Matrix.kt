package matrix

import java.lang.IllegalArgumentException
import java.lang.Math.pow
import kotlin.math.abs
import kotlin.math.sqrt

open class Matrix constructor(row: Int,
                              column: Int = row,
                              elements: Array<DoubleArray> = Array<DoubleArray>(row) { DoubleArray(column) }) {
    var row = row
        protected set
    var column = column
        protected set
    var elements = elements
        protected set
    val isSquare: Boolean get() = row == column

    //constructor(row:Int,column:Int,elements:Array<DoubleArray>):this(row,column){this.elements=elements}
    //constructor(n:Int):this(n,n)
    //constructor(n:Int,elements:Array<DoubleArray> = Array<DoubleArray>(n){DoubleArray(n)}):this(n,n,elements)
    operator fun plus(other: Matrix): Matrix {
        if (this.row == other.row && this.column == other.column) {
            val result = Matrix(row, column)
            for (i in 0 until row)
                for (j in 0 until column)
                    result.elements[i][j] = elements[i][j] + other.elements[i][j]
            return result
        } else throw IllegalArgumentException("The two matrices must be identical in addition and subtraction!")
    }

    operator fun minus(other: Matrix): Matrix {
        if (this.row == other.row && this.column == other.column) {
            val result = Matrix(row, column)
            for (i in 0 until row)
                for (j in 0 until column)
                    result.elements[i][j] = elements[i][j] - other.elements[i][j]
            return result
        } else throw IllegalArgumentException("The two matrices must be identical in addition and subtraction!")
    }

    //operator fun plusAssign(other:Double) = elements.forEach { it.forEach { x:Double->x+=other } }
    operator fun plus(other: Double): Matrix {
        val result = Matrix(row, column)
        for (i in 0 until row)
            for (j in 0 until column)
                result.elements[i][j] = elements[i][j] + other
        return result
    }

    operator fun plusAssign(other: Double) {
        for (i in 0 until row)
            for (j in 0 until column)
                elements[i][j] += other
    }

    operator fun plusAssign(other: Matrix) {
        if (this.row == other.row && this.column == other.column) {
            for (i in 0 until row)
                for (j in 0 until column)
                    elements[i][j] += other.elements[i][j]
        } else throw IllegalArgumentException("The two matrices must be identical in addition and subtraction!")
    }

    operator fun minusAssign(other: Matrix) {
        if (this.row == other.row && this.column == other.column) {
            for (i in 0 until row)
                for (j in 0 until column)
                    elements[i][j] -= other.elements[i][j]
        } else throw IllegalArgumentException("The two matrices must be identical in addition and subtraction!")
    }

    open fun transpose(): Matrix {
        val result = Matrix(column, row)
        for (i in 0 until row)
            for (j in 0 until column)
                result.elements[j][i] = elements[i][j]
        return result
    }

    operator fun times(other: Matrix): Matrix {
        if (column == other.row) {
            val result = Matrix(row, other.column)
            for (i in 0 until row)
                for (j in 0 until other.column)
                    for (k in 0 until column)
                        result.elements[i][j] += elements[i][k] * other.elements[k][j]
            return result
        } else throw IllegalArgumentException("The number of columns in the left matrix must equal to the number of rows in the right matrix!")
    }

    operator fun times(other: ColumnVector): ColumnVector {
        if (column == other.row) {
            val result = ColumnVector(DoubleArray(row))
            for (i in 0 until row)
                for (k in 0 until column)
                    result[i] += elements[i][k] * other[k]
            return result
        } else throw IllegalArgumentException("The number of columns in the left matrix must equal to the number of rows in the right matrix!")
    }

    operator fun timesAssign(other: Matrix): Unit {
        if (column == other.row) {
            val result: Array<DoubleArray> = Array<DoubleArray>(row) { DoubleArray(column) }
            column = other.column
            for (i in 0 until row)
                for (j in 0 until other.column)
                    for (k in 0 until column)
                        result[i][j] += elements[i][k] * other.elements[k][j]
            elements = result
        } else throw IllegalArgumentException("The number of columns in the left matrix must equal to the number of rows in the right matrix!")
    }

    operator fun timesAssign(other: Double): Unit {
        for (i in 0 until row)
            for (j in 0 until column)
                elements[i][j] *= other
    }

    open operator fun times(other: Double): Matrix {
        val result = Matrix(row, column)
        for (i in 0 until row)
            for (j in 0 until column)
                result.elements[i][j] = elements[i][j] * other
        return result
    }

    open operator fun div(other: Double): Matrix {
        val result = Matrix(row, column)
        for (i in 0 until row)
            for (j in 0 until column)
                result.elements[i][j] = elements[i][j] / other
        return result
    }

    operator fun unaryMinus(): Matrix {
        val result = Matrix(row, column)
        for (i in 0 until row)
            for (j in 0 until column)
                result.elements[i][j] = -elements[i][j]
        return result
    }
}

class SquareMatrix(val size: Int,
                   elements: Array<DoubleArray> = Array<DoubleArray>(size) { DoubleArray(size) }) : Matrix(size, size, elements)

interface Vector {
    val size: Int //get() = vector.size
    operator fun get(index: Int): Double //= vector[index]
    operator fun set(index: Int, value: Double) /*{
        vector[index] = value
    }*/
    fun norm(p: Any = 2): Double {
        var result = 0.0
        for (i in 0 until size)
            when (p) {
                0, 0.0 -> result += if (this[i] == 0.0) 0 else 1
                1 -> result += abs(this[i])
                2 -> result += this[i] * this[i]
                is Double -> result += pow(abs(this[i]), p)
                "Infinity" -> if (abs(this[i]) > result) result = abs(this[i])
            }
        return when (p) {
            0, 0.0, 1, "Infinity" -> result
            2 -> sqrt(result)
            is Double -> pow(result, 1 / p)
            else -> throw IllegalArgumentException("P must be number or Infinity")
        }
    }
}

class ColumnVector(vector: DoubleArray) : Matrix(vector.size, 1, Array(vector.size) { i -> doubleArrayOf(vector[i]) }), Vector {
    override val size: Int = row
    override operator fun get(index: Int) = elements[index][0]
    override operator fun set(index: Int, value: Double) {
        elements[index][0] = value
    }

    //constructor(m: Matrix, index: Int = 0) : this(DoubleArray(m.row) { i->m.elements[i][index] })

    operator fun plus(other: ColumnVector) = ColumnVector(DoubleArray(size) { i -> this[i] + other[i] })
    operator fun minus(other: ColumnVector) = ColumnVector(DoubleArray(size) { i -> this[i] - other[i] })
    override operator fun times(other: Double) = ColumnVector(DoubleArray(size) { i -> this[i] * other })
    override operator fun div(other: Double) = ColumnVector(DoubleArray(size) { i -> this[i] / other })
    override fun transpose(): RowVector = RowVector(DoubleArray(size) { i -> this[i] })
}


class RowVector(private val vector: DoubleArray) : Matrix(1, vector.size, arrayOf(vector)), Vector {
    override val size: Int get() = column
    override operator fun get(index: Int) = vector[index]
    override operator fun set(index: Int, value: Double) {
        vector[index] = value
    }

    override fun transpose(): ColumnVector = ColumnVector(vector.clone())
}

fun Diagonal(diag: DoubleArray): Matrix {
    val n = diag.size
    val result = SquareMatrix(n)
    for (i in 0 until n)
        result.elements[i][i] = diag[i]
    return result
}

fun Identity(n: Int): Matrix = Diagonal(DoubleArray(n) { 1.0 })

fun println(m: Matrix) = m.elements.forEach { it.forEach { print("$it ") };println() }

fun dot(A: Vector, B: Vector): Double {
    var result: Double = 0.0
    for (i in 0 until A.size)
        result += A[i] * B[i]
    return result
}


