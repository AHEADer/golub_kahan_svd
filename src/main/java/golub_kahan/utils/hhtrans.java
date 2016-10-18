package golub_kahan.utils;

/**
 * Created by david on 10/10/16.
 */
import golub_kahan.datatype.*;
import org.apache.log4j.Logger;
import scala.None;

import java.util.ArrayList;
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
        Matrix P_original = new Matrix(A.row, A.row);
        Matrix Q_original = new Matrix(A.col, A.col);
        List<Matrix> return_list = new ArrayList<>();
        List<Vec> vec_for_P = A.cutByColOrRow(true);
        List<Vec> vec_for_Q = A.cutByColOrRow(false);

        return return_list;
    }
}
