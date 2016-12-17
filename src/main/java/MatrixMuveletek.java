import static java.lang.Math.*;

public class MatrixMuveletek {

    private static final double EPSILON = 1e-10;

    private static void osszeadas(double A[][],double B[][]) {

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

    private static void kivonas(double A[][],double B[][]) {

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

    private static void konstansSzorzas(double A[],double B[][]) {

        double R[][] = new double[2][2];

        if( A.length == 1 && x2x2(B)) {

            for(int i = 0; i < 2; i++) {
                for(int j=0; j < 2; j++) {
                    R[i][j] = A[0] * B[i][j];
                }
            }
            printMatrix(R);
        } else {
            System.out.println("Nem megfelelo matrixok konstan x matrix muvelethez!");
        }
    }

    private static void mulMatrices(double A[][], double B[][]) {

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

    private static void inverse(double A[][]) {

        double subMatrix[][] = new double[2][2];

        if(x2x2(A)) {

            double multiplier[] = {1 / ( (A[0][0] * A[1][1]) - (A[0][1] * A[1][0]) )};

            subMatrix[0][0] = A[1][1];
            subMatrix[0][1] = - A[0][1];
            subMatrix[1][0] = - A[1][0];
            subMatrix[1][1] = A[0][0];

            konstansSzorzas(multiplier, subMatrix);

        } else {
            System.out.println("Nem 2x2 matrix.");
        }
    }

    private static void determinant(double A[][]) {

        double det;

        if(x2x2(A)) {

            det = ( (A[0][0] * A[1][1]) - (A[0][1] * A[1][0]) );

            System.out.println(det);
        } else {
            System.out.println("Nem 2x2 matrix.");
        }
    }

    private static void transponse(double A[][]) {

        double transponsed[][] = new double[A[0].length][A.length];

        if(x2x2(A)) {

            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A[0].length; j++) {
                    transponsed[j][i] = A[i][j];
                }
            }
            printMatrix(transponsed);
        } else {
            System.out.println("Nem 2x2 matrix.");
        }
    }

    private static void gaussElimination(double A[][], double b[]) {

        int N  = b.length;

        for (int p = 0; p < N; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (abs(A[i][p]) > abs(A[max][p])) {
                    max = i;
                }
            }
            double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            double   t    = b[p]; b[p] = b[max]; b[max] = t;

            // singular or nearly singular
            if (abs(A[p][p]) <= EPSILON) {
                throw new RuntimeException("Matrix is singular or nearly singular");
            }

            // pivot within A and b
            for (int i = p + 1; i < N; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < N; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // back substitution
        double[] x = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        // print results
        for (int i = 0; i < N; i++) {
            System.out.println(x[i]);
        }
    }

    private static void eigenvalues(double A[][]) {

        double eValues[] = new double[2];

        if( x2x2(A)) {

            double a = 1;
            double b = ( (-1) * (A[0][0]) ) - A[1][1];
            double c = (A[0][0] * A[1][1]) - (A[0][1] * A[1][0]);

            double temp1 = sqrt(b * b - 4 * a * c);

            eValues[0] = (-b +  temp1) / (2*a) ;
            eValues[1] = (-b -  temp1) / (2*a) ;
        }

        for(int i = 0; i < eValues.length; i++) {
            System.out.println(i + ". sajatertek " + eValues[i]);
        }
    }

    private static void eigenvectors(double A[][]){

        double eValues[] = new double[2];
        double eVectors[][] = new double[2][2];

        if( x2x2(A)) {

            double a = 1;
            double b = ( (-1) * (A[0][0]) ) - A[1][1];
            double c = (A[0][0] * A[1][1]) - (A[0][1] * A[1][0]);

            double temp1 = sqrt(b * b - 4 * a * c);

            eValues[0] = (-b +  temp1) / (2*a) ;
            eValues[1] = (-b -  temp1) / (2*a) ;
        }

        if(A[0][1] == 0 && A[1][0] == 0) {
            eVectors[0][0] = 1;
            eVectors[0][1] = 0;

            eVectors[1][0] = 0;
            eVectors[1][1] = 1;

            for (int i = 0; i < 2; i++) {
                System.out.println(i + 1 + ". eigenvector: ");
                for(int j = 0; j < 2; j++) {
                    System.out.println(eVectors[i][j]);
                }
                System.out.println();
            }

            return;
        }

        if(A[1][0] != 0) {
            eVectors[0][0] = eValues[0] - A[1][1];
            eVectors[0][1] = A[1][0];

            eVectors[1][0] = eValues[1] - A[1][1];
            eVectors[1][1] = A[1][0];

            for (int i = 0; i < 2; i++) {
                System.out.println(i + 1 + ". eigenvector: ");
                for(int j = 0; j < 2; j++) {
                    System.out.println(eVectors[i][j]);
                }
                System.out.println();
            }

            return;
        }

        if(A[0][1] != 0) {
            eVectors[0][0] = A[0][1];
            eVectors[0][1] = eValues[0] - A[0][0];

            eVectors[1][0] = A[0][1];
            eVectors[1][1] = eValues[1] - A[0][0];

            for (int i = 0; i < 2; i++) {
                System.out.println(i + 1 + ". eigenvector: ");
                for(int j = 0; j < 2; j++) {
                    System.out.println(eVectors[i][j]);
                }
                System.out.println();
            }

            return;
        }

    }

    private static void printMatrix(double A[][]) {
        for(int i = 0; i < A.length; i++) {
            for(int j=0; j < A[0].length; j++) {
                System.out.print(A[i][j] + ", ");
            }
            System.out.println();
        }
    }

    private static Boolean x2x2(double A[][]) {
        if( (A.length == 2 && A[0].length == 2) ){
            return true;
        } else {
            return false;
        }
    }

    private static Boolean canMultiply(double A[][], double B[][]) {
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

        double P5[][] = {
                {2, 0},
                {1, 1}
        };

        double K[] = {10};

        double gP1[][] = {
                {2, 7, 1},
                {4, 0, 3},
                {1, 3, 1}
        };

        double gP2[] = {1, -1, -1};

        double[][] gP3 = {
                { 1, 5, -2 },
                { 2, 3, 1 },
                { 2, 4, -3 }
        };

        double[] gP4 = { 2, 5, 2 };

        double eP[][] = {
                {0, 1},
                {-2, -3}
        };

        double eP2[][] = {
                {2, 2},
                {2, 2}
        };

        double eP3[][] = {
                {3, 0},
                {8, -1}
        };

        osszeadas(P1, P2);
        System.out.println();
        kivonas(P1, P2);
        System.out.println();
        konstansSzorzas(K, P1);
        System.out.println();
        mulMatrices(P1,P2);
        System.out.println();
        inverse(P5);
        System.out.println();
        determinant(P1);
        System.out.println();
        transponse(P5);
        System.out.println();
        gaussElimination(gP3, gP4);
        System.out.println();
        eigenvalues(eP2);
        System.out.println();
        eigenvectors(eP2);

    }
}