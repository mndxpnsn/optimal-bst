import java.util.Random;

public class BST {

    static int LARGE_NUM = 1000;

    static double[][][] dp;
    static int[][] root;

    static double min(double a, double b) {
        double res = 0;

        if(a < b) res = a;
        else res = b;

        return res;
    }

    static double optimal_bst(int i, int j, int d, double[] p, double[] q) {
        double res = LARGE_NUM;

        // Get results from memo table if available
        if(dp[i][j][d] != 0.0) {
            return dp[i][j][d];
        }

        // Case is a leaf
        if(i == j) {
            double val1 = (d + 1) * p[i];
            double val2 = (d + 2) * q[i];
            double val3 = (d + 2) * q[i + 1];
            root[i][j] = i;
            return val1 + val2 + val3;
        }

        // Case is subtree
        if(j > i) {
            double val1 = (d + 1) * p[i] + (d + 2) * q[i] + optimal_bst(i + 1, j, d + 1, p, q);
            double val2 = (d + 1) * p[j] + (d + 2) * q[j + 1] + optimal_bst(i, j - 1, d + 1, p, q);
            for(int k = i + 1; k <= j - 1; ++k) {
                double val1_loc = optimal_bst(i, k - 1, d + 1, p, q);
                double val2_loc = optimal_bst(k + 1, j, d + 1, p, q);
                double res_loc = (d + 1) * p[k] + val1_loc + val2_loc;

                // Get root
                if(res > res_loc) {
                    res = res_loc;
                    root[i][j] = k;
                }
            }

            // Case node i is (sub) root
            if(res > val1) {
                res = val1;
                root[i][j] = i;
            }

            // Case node j is (sub) root
            if(res > val2) {
                res = val2;
                root[i][j] = j;
            }
        }

        // Store results in memo table
        dp[i][j][d] = res;

        return res;
    }

    static double[][] optimal_bst_ref(double[] p, double[] q, int n) {
        double[][] e = new double[n + 2][n + 2];
        double[][] w = new double[n + 2][n + 2];

        for(int i = 1; i <= n + 1; ++i) {
            e[i][i - 1] = q[i - 1];
            w[i][i - 1] = q[i - 1];
        }

        for(int l = 1; l <= n; ++l) {
            for(int i = 1; i <= n - l + 1; ++i) {
                int j = i + l - 1;
                e[i][j] = LARGE_NUM;
                w[i][j] = w[i][j - 1] + p[j] + q[j];
                for(int r = i; r <= j; ++r) {
                    double t = e[i][r - 1] + e[r + 1][j] + w[i][j];
                    if(t < e[i][j]) {
                        e[i][j] = t;
                    }
                }
            }
        }

        return e;
    }

    static void init_vec(int n, double[] v) {
        // Initialize probability vector v with random data.
        // Note the probabilities do not sum to one.
        for(int i = 0; i < n; ++i) {
            Random r = new Random();
            double min = 0.0;
            double max = 1.0 / (2.0 * n);

            double rand_val = min + (max - min) * r.nextDouble();

            v[i] = rand_val;
        }
    }

    static void set_vec(double[] vi, int n, double[] v) {

        v[0] = 0.0;

        for(int i = 1; i < n + 1; ++i) {
            v[i] = vi[i - 1];
        }
    }

    static void print_optimal_bst(int[][] root, int i, int j, int n) {

        // Compute root key
        int k = root[i][j];

        // Print root key
        if(j - i == n - 1) {
            System.out.println("root key is: " + "k" + (k + 1));
        }

        // Check if at left boundary
        if(k == i) {
            System.out.println("d" + (i) + " is the left child of k" + (k + 1));
        }

        // Check if at right boundary
        if(k == j) {
            System.out.println("d" + (j + 1) + " is the right child of k" + (k + 1));
        }

        // Recursively print optimal bst
        if(j > i) {

            if(k - 1 >= i) {
                int l = root[i][k - 1];
                System.out.println("k" + (l + 1) + " is the left child of k" + (k + 1));
                print_optimal_bst(root, i, k - 1, n);
            }

            if(k + 1 <= j) {
                int r = root[k + 1][j];
                System.out.println("k" + (r + 1) + " is the right child of k" + (k + 1));
                print_optimal_bst(root, k + 1, j, n);
            }
        }
    }

    public static void main(String[] args) {

        // Number of nodes in bst
        int n = 15;

        // Declare and initialize probability vectors with random data
        double[] p = new double[n];
        double[] q = new double[n + 1];
        double[] p_ref = new double[n + 1];

        init_vec(n, p);
        init_vec(n + 1, q);
        set_vec(p, n, p_ref);

        // Initialize memo and root tables
        dp = new double[n][n][n];
        root = new int[n + 1][n + 1];

        // Compute optimal bst
        double opt_bst = optimal_bst(0, n - 1, 0, p, q);

        // Compute optimal bst using reference method
        double[][] e = optimal_bst_ref(p_ref, q, n);

        // Print results
        print_optimal_bst(root, 0, n - 1, n);

        System.out.println("optimal bst cost: " + opt_bst);
        System.out.println("optimal bst cost reference method: " + e[1][n]);
    }
}