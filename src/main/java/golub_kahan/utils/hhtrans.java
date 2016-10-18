package golub_kahan.utils;

/**
 * Created by david on 10/10/16.
 */
import golub_kahan.datatype.*;
import org.apache.log4j.Logger;
import scala.None;

import java.util.List;

/*household transformation*/
public class hhtrans {
    private final Logger log;
    private int tmp_a;

    public hhtrans(Logger log) {
        this.log = log;
    }

/*
A is a m*n matrix, and it can be decomposed as A=PJQ.
P and Q are unitary matrices and J is an m*n bidiagonal.
returnlist cantains three matrix PJQ. P is a m*m matrix
and Q is a n*n matrix.
*/
    public List<Matrix> HouseholdTransformation(Matrix A) {
        Matrix P_original;
        return None;
    }
}
