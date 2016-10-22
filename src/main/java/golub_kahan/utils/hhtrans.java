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
    public static List<Matrix> HouseholdTransformation(Matrix A) {
        Matrix P_original = new Matrix(A.row, A.row);
        Matrix Q_original = new Matrix(A.col, A.col);
        List<Matrix> return_list = new ArrayList<>();
        List<Vec> vec_for_P = A.cutByColOrRow(true);
        List<Vec> vec_for_Q = A.cutByColOrRow(false);
        A.show();
        for (int i = 0;i < vec_for_P.size();i++) {
            Vec currentVec = vec_for_P.get(i);
            int size = currentVec.getSize();
            Vec UnitVec = currentVec.formUnitVec(i);
            double sigma = currentVec.norm();
            Vec V = vecAdd(currentVec, UnitVec.multiply(sigma));
            V = V.divide(V.norm());
            Matrix H = new Matrix(size, size);
            H.formIdentityMatrix();
            H = matrixSubtract(H, matrixMultiply(2,outProduct(V, V)));

        }
        return return_list;
    }

    public static Vec vecAdd(Vec a, Vec b) {
        if (a.getSize()!=b.getSize())
            throw new IllegalArgumentException();
        Vec result = new Vec(a.getSize());
        for (int j = 0; j < a.getSize(); j++) {
            result.set(j, a.get(j)+b.get(j));
        }
        return result;
    }

    public static Matrix matrixSubtract(Matrix a, Matrix b) {
        if (a.col!=b.col||a.row!=b.row)
            throw new IllegalArgumentException();
        Matrix result = new Matrix(a.row, a.col);
        for (int i = 0; i < a.row; i++) {
            for (int j = 0; j < a.col; j++) {
                result.set(i, j, a.get(i, j) - b.get(i, j));
            }
        }
        return result;
    }

    public static Matrix matrixAddition(Matrix a, Matrix b) {
        if (a.col!=b.col||a.row!=b.row)
            throw new IllegalArgumentException();
        Matrix result = new Matrix(a.row, a.col);
        for (int i = 0; i < a.row; i++) {
            for (int j = 0; j < a.col; j++) {
                result.set(i, j, a.get(i, j) + b.get(i, j));
            }
        }
        return result;
    }

    public static Matrix matrixMultiply(double a, Matrix b) {
        Matrix result = new Matrix(b.row, b.col);
        for (int i = 0; i < b.row; i++) {
            for (int j = 0; j < b.col; j++) {
                result.set(i, j, a * b.get(i, j));
            }
        }
        return result;
    }

    //a*bT
    public static Matrix outProduct(Vec a, Vec b) {
        Matrix result = new Matrix(a.getSize(), b.getSize());
        for (int i = 0; i < a.getSize(); i++) {
            for (int j = 0; j < b.getSize(); j++) {
                result.set(i, j, a.get(i) * b.get(j));
            }
        }
        return result;
    }
    public static Vec renderVec(Matrix a, Vec b) {
        if (a.row!=)
    }
}
