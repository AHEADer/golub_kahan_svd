package golub_kahan;

import golub_kahan.datatype.Matrix;
import golub_kahan.utils.hhtrans;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

public class entrance {

    public static void main(String[] args) throws InterruptedException {
        Matrix A = new Matrix(10, 10);
        A.random(10, 10);
        Logger log = Logger.getLogger(hhtrans.class);
        hhtrans test = new hhtrans(log);
        List<Matrix> PJQ = new ArrayList<>();
        PJQ = hhtrans.HouseholdTransformation(A);
        PJQ.get(1).show();
    }
}
