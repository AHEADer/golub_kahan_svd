package golub_kahan.datatype;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;

import java.io.*;
import java.util.*;

/**
 * Created by cxy on 15-10-21.
 */
public class Matrix implements Serializable{
    private double matrix[][];
    public int row;
    public int col;



    public Matrix(int row, int col) {
        this.row = row;
        this.col = col;
        matrix = new double[row][col];
    }

    public double get(int i, int j) { return matrix[i][j]; }

    public double[][] getMatrix() {
        return matrix;
    }

    public double norm() {
        double result = 0;
        for(double[] arr : matrix) {
            for(double value : arr) {
                result += value*value;
            }
        }
        return Math.sqrt(result);
    }

    public boolean isEmpty() {
        if(row == 0 && col == 0)
            return true;
        else
            return false;
    }

    public void set(int i, int j, double value) {
        matrix[i][j] = value;
    }

    //根据一个列向量来生成一个矩阵的一列
    public void set(int colIndex, Vec arr) {
        for(int i = 0; i < arr.getSize(); i++) {
            matrix[i][colIndex] = arr.get(i);
        }
    }

    public Vec getCol(int colIndex) {
        ArrayList<Double> result = new ArrayList<>();
        for(int i = 0; i < row; i++) {
            result.add(matrix[i][colIndex]);
        }
        return Vec.fromList(result);
    }

    public Vec getRow(int rowIndex) {
        return Vec.fromArray(matrix[rowIndex]);
    }


    public Matrix trans() {
        Matrix result = new Matrix(col,row);
        for(int i = 0; i < col; i++)
            for(int j = 0; j < row; j++)
                result.set(i,j,get(j,i));

        return result;
    }

    public Matrix dot(Matrix that) {
        if(this.col != that.row) {
            System.out.println("Matrix dot operation is not matched!");
            throw new IllegalArgumentException();
        } else {
            Matrix result = new Matrix(this.row,that.col);
            for(int i = 0; i < result.row; i++) {
                for(int j = 0; j < result.col; j++) {
                    Vec rowVec = getRow(i);
                    Vec colVec = that.getCol(j);
                    result.set(i,j,rowVec.dot(colVec));
                }
            }
            return result;
        }
    }

    //矩阵乘以一个对角矩阵，对角矩阵用数组表示
    //用以节省存储空间
    public void dotInSelf(double[] SIGMA) {
        if(this.col != SIGMA.length) {
            System.out.println("Matrix dot operation is not matched!");
            throw new IllegalArgumentException();
        } else {
            for(int i = 0; i < this.col; i++) {
                for(int j = 0; j < this.row; j++) {
                    this.set(j,i, get(j,i)*SIGMA[i]);
                }
            }
        }
    }

    //矩阵乘以一个对角矩阵，对角矩阵用数组表示
    //用以节省存储空间
    public Matrix dot(double[] SIGMA) {
        if(this.col != SIGMA.length) {
            System.out.println("Matrix dot operation is not matched!");
            throw new IllegalArgumentException();
        } else {
            Matrix result = new Matrix(this.row,this.col);
            for(int i = 0; i < result.col; i++) {
                for(int j = 0; j < result.row; j++) {
                    result.set(j,i, get(j,i)*SIGMA[i]);
                }
            }
            return result;
        }
    }

    public Matrix minus(Matrix that) {
        if( (this.row != that.row) && (this.col != that.col) ) {
            System.out.println("Matrix dot operation is not matched!");
            throw new IllegalArgumentException();
        } else {
            Matrix result = new Matrix(this.row,this.col);
            for(int i = 0; i < result.row; i++) {
                for(int j = 0; j < result.col; j++) {
                    result.set(i,j,Math.abs(Math.abs(this.get(i, j)) - Math.abs(that.get(i, j))));
                }
            }
            return result;
        }
    }

    @Override
    public String toString() {
//        String result = "[\n";
        String result = "";
        for(double[] arr: matrix) {
            for(double value : arr){
                result += value;
                result += "\b\t";
            }
            result += "\n";
        }
//        result += "]\n";
        return result;
    }

    public static class TupleComparator implements Comparator<Tuple> {

        @Override
        public int compare(Tuple o1, Tuple o2) {
            if (o1 == o2) {
                return 0;
            }
            else {
                if (o1.getNorm() > o2.getNorm())
                    return -1;
                else if (o1.getNorm() < o2.getNorm())
                    return 1;
                else
                    return 0;
            }
        }
    }

    public static class Tuple {
        private double norm;
        private Vec vector;
        public Tuple(double norm, Vec vector) {
            this.norm = norm;
            this.vector = vector;
        }

        public double getNorm(){ return norm; }
        public Vec getVector() { return vector; }
    }

    public static Matrix random(int row, int col) {
        Random rand = new Random(47);
        Matrix result = new Matrix(row,col);

        for(int i = 0; i < row; i++)
            for(int j = 0; j < col; j++)
                result.set(i,j,rand.nextDouble());

        return result;
    }

    public static Matrix from2DArray(double[][] arr) {
        int row = arr.length;
        int col = arr[0].length;
        Matrix result = new Matrix(row,col);

        for(int i = 0; i < row; i++)
            for(int j = 0; j < col; j++)
                result.set(i,j,arr[i][j]);

        return result;
    }

    public static Matrix emptyMatrix() {
        return new Matrix(0,0);
    }

    public static Matrix eye(int size) {
        Matrix result = new Matrix(size,size);

        for(int i = 0; i < size; i++) {
            result.set(i,i,1);
        }

        return result;
    }

    public static Matrix readMatrix(String path) throws IOException{
        File file;
        CSVReader csvReader = null;
        Matrix result = null;

        try {
            file = new File(path+".csv");
            csvReader = new CSVReader(new FileReader(file));
            List<String[]> data = csvReader.readAll();

            int x = data.size();
            int y = data.get(0).length;

            result = new Matrix(x,y);

            for (int i = 0; i < x; i++) {
                for (int j = 0; j < y; j++) {
                    result.set(i, j, Double.parseDouble(data.get(i)[j]));
//                    System.out.print(Double.parseDouble(data.get(i)[j]) + " ");
                }
            }
        } finally {
            csvReader.close();
        }

        return result;
    }

    //将矩阵写入到文件中，文件名为filename
    public  static void writeMatrixToFile(Matrix matrix, String filename) throws IOException {
        File path = new File(filename+".csv");
        Writer writer = new FileWriter(path);
        CSVWriter csvWriter = new CSVWriter(writer);

        String[] temp = new String[matrix.getMatrix()[0].length];
        for(double[] x : matrix.getMatrix()) {
            for (int i = 0; i < temp.length; i++) {
                temp[i] = String.valueOf(x[i]);
            }
            csvWriter.writeNext(temp);
        }
        csvWriter.close();
    }

    //if mode is True,combine by col; if False, combine by row
    public void combineByColOrRow(List<Vec> vec_col_list, boolean mode) {
        int judge = mode?col:row;
        for (int i = 0;i < vec_col_list.size();i++) {
            if (vec_col_list.get(i).getSize()!=judge) {
                System.out.println("The "+i+"th Vec is not compatible with the Matrix");
                throw new IllegalArgumentException();
            }
            else if (mode){
                for (int j = 0;j < judge;j++) {
                    matrix[j][i] = vec_col_list.get(i).get(j);
                }
            }
            else {
                for (int j = 0;j < judge;j++) {
                    matrix[i][j] = vec_col_list.get(i).get(j);
                }
            }
        }
    }

    public void show() {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                System.out.println(matrix[i][j]);
            }
        }
    }

    //if mode is True,cut by col; if False, cut by row
    public List<Vec> cutByColOrRow(boolean mode) {
        List<Vec> vec_list = new ArrayList<>();
        if (mode) {
            for (int i = 0; i < col; i++) {
                Vec tmp = getCol(i);
                vec_list.add(i, tmp);
            }
        }
        else {
            for (int j = 0; j < row; j++) {
                Vec tmp2 = getRow(j);
                vec_list.add(j, tmp2);
            }
        }
        return vec_list;
    }

    public void formIdentityMatrix() {
        if (col!=row)
            throw new IllegalArgumentException();
        else {
            for (int i = 0; i < row; i++)
                for (int j = 0; j < col; j++) {
                    if (i == j)
                        matrix[i][j] = 1;
                    else
                        matrix[i][j] = 0;

                }
        }
    }

    public Matrix constructHermitianMatrix(List<Vec> refer_vec_list, boolean mode) {
        Matrix result = new Matrix(1,1);
        return result;
    }

    public static void main(String[] args) throws IOException, InterruptedException {

    }

}
