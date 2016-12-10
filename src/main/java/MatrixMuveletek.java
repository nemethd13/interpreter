

public class MatrixMuveletek {

    public static void osszeadas(double A[][],double B[][]) {

        double R[][] = new double[2][2];

        if( (A.length == 2 && A[0].length == 2) && (B.length == 2 && B[0].length == 2)) {

            for(int i = 0; i < 2; i++) {
                for(int j=0; j < 2; j++) {
                    R[i][j] = A[i][j] + B[i][j];
                }
            }
        } else {
            System.out.println("Nem 2x2-es matrixok!");
            return;
        }
        printMatrix(R);
    }

    public static void kivonas(double A[][],double B[][]) {

        double R[][] = new double[2][2];

        if( x2x2(A) && x2x2(B))  {

            for(int i = 0; i < 2; i++) {
                for(int j=0; j < 2; j++) {
                    R[i][j] = A[i][j] - B[i][j];
                }
            }
        } else {
            System.out.println("Nem 2x2-es matrixok!");
            return;
        }
        printMatrix(R);
    }

    public static void konstansSzorzas(double A[],double B[][]) {

        double R[][] = new double[2][2];

        if( A.length == 1 && x2x2(B)) {

            for(int i = 0; i < 2; i++) {
                for(int j=0; j < 2; j++) {
                    R[i][j] = A[0] * B[i][j];
                }
            }
        } else {
            System.out.println("Nem megfelelo matrixok konstan x matrix muvelethez!");
            return;
        }
        printMatrix(R);
    }

    public static void mulMatrices(double A[][], double B[][]) {

        double sum;
        double R[][] = new double[A.length][B[0].length];

        if( canMultiply(A, B)){
            for(int i=0; i < A.length; i++){
                for(int j=0; j < B[0].length; j++) {
                    sum = 0;
                    for(int k = 0; k < B.length; k++) {
                        sum += A[i][k] * B[k][j];
                    }
                    R[i][j] = sum;
                }
            }
            printMatrix(R);
        } else {
            System.out.println("A ket matrix nem osszeszorozhato!");

        }
    }

    public static void inverse(double A[][]) {

        if(x2x2(A)) {
            double konstans = 1 / ( (A[0][0] * A[1][1]) - 
        }
    }

    public static void printMatrix(double A[][]) {
        for(int i = 0; i < A.length; i++) {
            for(int j=0; j < A[0].length; j++) {
                System.out.print(A[i][j] + ", ");
            }
            System.out.println();
        }
    }

    public static Boolean x2x2(double A[][]) {
        if( (A.length == 2 && A[0].length == 2) ){
            return true;
        } else {
            return false;
        }
    }

    public static Boolean canMultiply(double A[][], double B[][]) {
        if( A[0].length == B.length) {
            return true;
        } else {
            return false ;
        }
    }

    public static void main(String[] args) {

        double P1[][] = {
                        {1, 2},
                        {2, 1}
        };

        double P2[][] = {
                        {2, 1},
                        {1, 2}
        };

        double P3[][] = {
                        {1, 2, 4},
                        {2, 1, 4}
        };

        double P4[][] = {
                        {1, 2, 3},
                        {1, 2, 3},
                        {1, 2, 3}
        };

        double K[] = {10};

        osszeadas(P1, P2);
        System.out.println();
        kivonas(P1, P2);
        System.out.println();
        konstansSzorzas(K, P1);
        System.out.println();
        mulMatrices(P1,P2);

    }
}