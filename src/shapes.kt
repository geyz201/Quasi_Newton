package geometry.shapes
import Rectangle
import java.util.Random

fun createRandomRectangle():Rectangle{
    val random=Random();
    return Rectangle(random.nextInt(),random.nextInt())
}