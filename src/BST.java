import java.util.Random;

public class BST {

    static int LARGE_NUM = 1000;

    static double[][][] dp;

    static class OBST {
        static double val;
        static int[][] root;
    }

    static class OBST2 {
        static double val;
        static int[][] root;
    }

    static double min(double a, double b) {
        double res = 0;

        if(a < b) res = a;
        else res = b;

        return res;
    }

    static double opt_bst(int i, int j, int d, double[] p, double[] q, int[][] root) {
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
            double val1 = (d + 1) * p[i] + (d + 2) * q[i] + opt_bst(i + 1, j, d + 1, p, q, root);
            double val2 = (d + 1) * p[j] + (d + 2) * q[j + 1] + opt_bst(i, j - 1, d + 1, p, q, root);
            for(int k = i + 1; k <= j - 1; ++k) {
                double val1_loc = opt_bst(i, k - 1, d + 1, p, q, root);
                double val2_loc = opt_bst(k + 1, j, d + 1, p, q, root);
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

    static OBST2 optimal_bst_ref(double[] p, double[] q, int n) {
        OBST2 obst = new OBST2();
        double[][] e = new double[n + 2][n + 2];
        double[][] w = new double[n + 2][n + 2];
        obst.root = new int[n + 2][n + 2];

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
                        obst.root[i][j] = r;
                    }
                }
            }
        }

        obst.val = e[1][n];

        return obst;
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

    static void print_optimal_bst_ref(int[][] root, int i, int j, int n) {

        // Compute root key
        int k = root[i][j];

        // Print root key
        if(j - i == n - 1) {
            System.out.println("root key is: " + "k" + (k));
        }

        // Check if at left boundary
        if(k == i) {
            System.out.println("d" + (i - 1) + " is the left child of k" + (k));
        }

        // Check if at right boundary
        if(k == j) {
            System.out.println("d" + (j) + " is the right child of k" + (k));
        }

        // Recursively print optimal bst
        if(j > i) {

            if(k - 1 >= i) {
                int l = root[i][k - 1];
                System.out.println("k" + (l) + " is the left child of k" + (k));
                print_optimal_bst_ref(root, i, k - 1, n);
            }

            if(k + 1 <= j) {
                int r = root[k + 1][j];
                System.out.println("k" + (r) + " is the right child of k" + (k));
                print_optimal_bst_ref(root, k + 1, j, n);
            }
        }
    }

    static OBST optimal_bst(double[] p, double[] q, int n) {

        // Initialize memo and root tables
        OBST obst = new OBST();
        obst.root = new int[n][n];
        dp = new double[n][n][n];

        // Compute optimal binary search tree
        obst.val = opt_bst(0, n - 1, 0, p, q, obst.root);

        return obst;
    }

    public static void main(String[] args) {

        // Number of nodes in bst
        int n = 10;

        // Declare and initialize probability vectors with random data
        double[] p = new double[n];
        double[] q = new double[n + 1];
        double[] p_ref = new double[n + 1];

        init_vec(n, p);
        init_vec(n + 1, q);
        set_vec(p, n, p_ref);

        // Compute optimal bst
        OBST obst = optimal_bst(p, q, n);

        // Compute optimal bst using reference method
        OBST2 obst_ref = optimal_bst_ref(p_ref, q, n);

        // Print results
        System.out.println("Printing optimum bst");
        print_optimal_bst(obst.root, 0, n - 1, n);

        System.out.println("Printing optimum bst reference method");
        print_optimal_bst_ref(obst_ref.root, 1, n, n);

        System.out.println("optimal bst cost: " + obst.val);
        System.out.println("optimal bst cost reference method: " + obst_ref.val);
    }
}