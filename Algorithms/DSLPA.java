import java.io.*;
import java.util.*;

class Node {
    String id = "";//Record node number
    ArrayList<Node> neighbor = new ArrayList<Node>();//Store neighbor nodes
    ArrayList<Node> neighborAndMe = new ArrayList<Node>();//Storage node's neighbor nodes and itself
    int tag = 0;// label
    int n = 0;// degree
    double p = 0;// Local density
    double d = 0;// distance
    double px = 0;// standard deviation normalized local density
    double dx = 0;// standard deviation normalized distance
    double y = 0;// y=px*dx
    public double[] score = new double[]{1.0, 0};//Record node importance
    double inf = 0;//Node influence
}

class Edge {
    String source = "";
    String target = "";
}

public class DSLPA {
    public ArrayList<Node> node = new ArrayList<>(); // Storage node set
    public ArrayList<Edge> edge = new ArrayList<>();// Storage edge set
    public double[][] jaccard; // Jaccard coefficient matrix
    public double epsilon = 2.0;//The ε parameter in DS-LPA needs to be set manually
    public int cid = 1;// Record the number of the community
    public double infmax = 0;//Record the maximum influence in the network
    public double infmin = 10000;//Record the minimum value of influence in the network
    public ArrayList<Node> center = new ArrayList<Node>();// Storage center node
    public int[] result; // Character-based coding method records experimental community partition(result[10]=1 means The community number of the node with id 10 is 1)
    public  int[] com; // Character-based coding method records common community partition
    public  int c_num = 0; // Number of communities for common community partition

    public static void main(String[] args) {
        DSLPA algorithm = new DSLPA();
        algorithm.dsLpa();
    }

    public void dsLpa() {
        read_real_dataset(node, edge, "karate.gml");//Read data set. Take Karate network as an example
        long startTime = System.currentTimeMillis();
        findCenter();
        PageRank();
        laberPropagating();
        long endTime = System.currentTimeMillis();
        try {
            readComFile("karate_comm.dat");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        //Print results
        System.out.println("NMI :"+NMI());
        System.out.println("ARI :"+ARI());
        System.out.println("Modularity :"+Q());
        System.out.println("Program running time:" + (endTime - startTime) + "ms");
    }

    /*
     * Read data set
     */
    public void read_real_dataset(ArrayList<Node> node, ArrayList<Edge> edge, String File) { 
        File file = new File(File);
        String L = "";
        int which = 0;
        String left = "";
        String right = "";
        Node nn = new Node();
        Edge ee = new Edge();
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            String s = null;
            while ((s = br.readLine()) != null) {
                L = s.trim();
                if (L.length() >= 4 && L.charAt(0) == 'n' && L.charAt(1) == 'o' && L.charAt(2) == 'd'
                        && L.charAt(3) == 'e') {
                    which = 0;
                    Node N = new Node();
                    nn = N;
                    node.add(nn);
                }
                if (L.length() >= 4 && L.charAt(0) == 'e' && L.charAt(1) == 'd' && L.charAt(2) == 'g'
                        && L.charAt(3) == 'e') {
                    which = 1;
                    Edge E = new Edge();
                    ee = E;
                    edge.add(ee);
                }
                if (which == 0) {
                    if (match(L, "["))
                        continue;
                    if (match(L, "]"))
                        continue;
                    if (L.length() == 0)
                        continue;
                    left = cutl(L);
                    right = cutr(L);
                    if (match(left, "id")) {
                        nn.id = right;
                    }
                }
                if (which == 1) {
                    if (match(L, "["))
                        continue;
                    if (match(L, "]"))
                        continue;
                    if (L.length() == 0)
                        continue;
                    left = cutl(L);
                    right = cutr(L);
                    if (match(left, "source")) {
                        ee.source = right;
                    }
                    if (match(left, "target")) {
                        ee.target = right;
                    }
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        if (node.get(0).id.equals("1")) {
            for (int i = 0; i < node.size(); i++) {
                node.get(i).id = Integer.parseInt(node.get(i).id) - 1 + "";
            }
            for (int i = 0; i < edge.size(); i++) {
                edge.get(i).source = Integer.parseInt(edge.get(i).source) - 1 + "";
                edge.get(i).target = Integer.parseInt(edge.get(i).target) - 1 + "";
            }
        }
        for (int i = 0; i < edge.size(); i++) {
            node.get((Integer.parseInt(edge.get(i).source))).neighbor
                    .add(node.get((Integer.parseInt(edge.get(i).target))));
            node.get((Integer.parseInt(edge.get(i).target))).neighbor
                    .add(node.get((Integer.parseInt(edge.get(i).source))));
        }
        for (int i = 0; i < node.size(); i++) {
            node.get(i).n = node.get(i).neighbor.size();
            node.get(i).neighborAndMe.add(node.get(i));
            node.get(i).neighborAndMe.addAll(node.get(i).neighbor);

        }
    }

    /*
     * String match
     */
    public boolean match(String la, String lb) {
        if (lb.length() != la.length())
            return false;
        for (int i = 0; i < lb.length(); i++) {
            if (la.charAt(i) == lb.charAt(i))
                continue;
            else
                return false;
        }
        return true;
    }

    /*
     * String processing:process target104 as target
     */
    public String cutl(String l) {
        String p = "";
        for (int i = 0; i < l.length(); i++) {
            if (l.charAt(i) != ' ')
                p = p + l.charAt(i);
            else
                break;
        }
        return p;
    }

    /*
     * String processing:process target104 as 104
     */
    public String cutr(String l) {
        String p = "";
        for (int i = 0; i < l.length(); i++) {
            if (l.charAt(i) != ' ')
                continue;
            else
                for (int j = i; j < l.length(); j++) {
                    p = p + l.charAt(j);
                }
        }
        p = p.trim();
        return p;
    }

    public void PageRank() {
        double bound = 0.0001;
        int oldIndex = 0;
        int newIndex = 1;
        int totaldegree = 0;
        for (int i = 0; i < node.size(); i++) {
            totaldegree += node.get(i).n;
        }
        int nodeNum = node.size();
        Node g = new Node();
        g.id = "-1";
        g.n = nodeNum;
        g.score[0] = 0;
        g.score[1] = 0;
        node.add(g);
        for (int iterNum = 0; iterNum < 100000; iterNum++) {
            for (int i = 0; i < node.size(); i++)
                node.get(i).score[newIndex] = 0;
            for (int i = 0; i < node.size(); i++) {
                if (i != node.size() - 1) {
                    for (int j = 0; j < node.get(i).neighbor.size(); j++)
                        node.get(i).score[newIndex] += (double) node.get(i).neighbor.get(j).score[oldIndex] / node.get(i).neighbor.get(j).n;
                    g.score[newIndex] += (double) node.get(i).score[oldIndex] / node.get(i).n;
                }
            }
            double diff = 0.0;
            for (int i = 0; i < node.size(); i++)
                diff += Math.abs(node.get(i).score[newIndex] - node.get(i).score[oldIndex]);
            oldIndex = newIndex;
            newIndex = 1 - oldIndex;
            if (diff < bound)
                break;
        }
        double avg = g.score[oldIndex] / nodeNum;
        for (int i = 0; i < node.size(); i++) {
            node.get(i).inf = node.get(i).score[oldIndex] + avg;
            if (node.get(i).inf > infmax)
                infmax = node.get(i).inf;
            if (node.get(i).inf < infmin)
                infmin = node.get(i).inf;
        }
        node.remove(g);
    }

    public void findCenter() {
        jaccard = new double[node.size()][node.size()];
        ArrayList<Node> copyList = new ArrayList<Node>(); // Temporary list for computing intersection
        double minp = 100000;// Record the minimum local density
        double maxp = 0;// Record the maximum value of local density
        double averagep = 0;// Mean of local density
        double mind = 100;// Minimum distance
        double maxd = 0;// Maximum distance
        double tempd;// Temporary variable used to calculate distance
        boolean bigFlag;// Mark for judging the point with the highest density
        double maxJaccard = 0;// Maximum Jaccard coefficient
        double sum = 0;// Used to calculate the sum of y
        double e = 0;// expectation of y
        double sd = 0;// standard deviation of y
        double temps = 0;// An intermediate quantity used to calculate the sum of standard deviations
        double bound;//Chebyshev upper bound

        // Calculate the Jaccard coefficient
        for (int i = 0; i < node.size(); i++) {
            for (int j = 0; j < node.size(); j++) {
                copyList = (ArrayList<Node>) node.get(i).neighborAndMe.clone();
                copyList.retainAll(node.get(j).neighborAndMe);
                jaccard[i][j] = (double) copyList.size()
                        / (node.get(i).neighborAndMe.size() + node.get(j).neighborAndMe.size() - copyList.size());
            }
        }

        // Take 1-Jacard coefficient as the distance between nodes
        for (int i = 0; i < node.size(); i++)
            for (int j = 0; j < node.size(); j++) {
                if (i == j) {
                    jaccard[i][j] = 0;
                    continue;
                }
                jaccard[i][j] = 1 - jaccard[i][j];
            }
        // Calculate the local density of each node
        for (int i = 0; i < node.size(); i++) {
            int s = 0;
            for (int j = 0; j < node.get(i).neighbor.size(); j++) {
                s = s + node.get(i).neighbor.get(j).n;
            }
            node.get(i).p = node.get(i).n + node.get(i).n * s * 1.0 / (node.get(i).n + s);
            if (node.get(i).p > maxp)
                maxp = node.get(i).p;
            if (node.get(i).p < minp)
                minp = node.get(i).p;
        }

        for (int i = 0; i < node.size(); i++) {
            node.get(i).px = (node.get(i).p - minp) / (maxp - minp);
            averagep += node.get(i).px;
        }
        averagep /= node.size();

        // Calculate the distance and its standard deviation
        for (int i = 0; i < node.size(); i++) {
            tempd = 100000;
            maxJaccard = 0;
            bigFlag = false;
            for (int j = 0; j < node.size(); j++) {
                if (j == i)
                    continue;
                if (jaccard[i][j] > maxJaccard)
                    maxJaccard = jaccard[i][j];
                if (node.get(j).px > node.get(i).px) {
                    bigFlag = true;
                    if (jaccard[i][j] < tempd)
                        tempd = jaccard[i][j];
                }
            }
            {
                if (bigFlag)
                    node.get(i).d = tempd;
                else
                    node.get(i).d = maxJaccard;
            }
            if (node.get(i).d > maxd)
                maxd = node.get(i).d;
            if (node.get(i).d < mind)
                mind = node.get(i).d;
        }
        for (int i = 0; i < node.size(); i++) {
            node.get(i).dx = (node.get(i).d - mind) / (maxd - mind);
        }

        /*
         * Calculate the constructor y, use Chebyshev's inequality to find the community center
         */
        for (int i = 0; i < node.size(); i++) {
            node.get(i).y = node.get(i).px * node.get(i).dx;
            sum = sum + node.get(i).y;
        }
        e = sum / node.size();
        for (int i = 0; i < node.size(); i++) {
            temps = temps + (node.get(i).y - e) * (node.get(i).y - e);
        }
        sd = Math.sqrt(temps / (node.size() - 1));
        bound = e + epsilon * sd;
        for (int i = 0; i < node.size(); i++) {
            if (node.get(i).y > bound) {
                center.add(node.get(i));
                node.get(i).tag = cid;
                cid++;
            }
        }

       /* Community initialization: Propagate the label of the community center to its neighbor nodes, the neighbor nodes need to meet
            1. Not a community center 2. Not connected to other community centers
       */
        for (int i = 0; i < center.size(); i++) {
            for (int j = 0; j < center.get(i).neighbor.size(); j++) {
                boolean centerflag = true;
                if (center.get(i).neighbor.get(j).tag == 0) {
                    for (int k = 0; k < center.get(i).neighbor.get(j).neighbor.size(); k++) {
                        if (center.contains(center.get(i).neighbor.get(j).neighbor.get(k))
                                && center.get(i).neighbor.get(j).neighbor.get(k).tag != center.get(i).tag)
                            centerflag = false;
                    }
                    if (centerflag)
                        center.get(i).neighbor.get(j).tag = center.get(i).tag;
                }
            }
        }

    }

    public void laberPropagating() {
        ArrayList<Node> update = new ArrayList<Node>();//Store the node whose label needs to be updated
        int neighborTagNum;//Record the number of labels in the neighbor
        double labelRate;
        double maxLabelRate;
        int updateIndex = 0;//Record the number of the node that needs to be updated

        // Store unlabeled nodes in update
        for (int i = 0; i < node.size(); i++) {
            if (node.get(i).tag == 0)
                update.add(node.get(i));
        }
        //    System.out.println(update.size());
        while (!update.isEmpty()) {
            maxLabelRate = 0;
            for (int i = 0; i < update.size(); i++) {
                neighborTagNum = 0;
                for (int j = 0; j < update.get(i).neighbor.size(); j++) {
                    if (update.get(i).neighbor.get(j).tag != 0)
                        neighborTagNum++;
                }
                labelRate = ((double) neighborTagNum / update.get(i).neighbor.size()) * (1 + (update.get(i).inf - infmin) / (infmax - infmin));
                if (labelRate > maxLabelRate) {
                    maxLabelRate = labelRate;
                    updateIndex = i;
                }
            }
            if (maxLabelRate == 0) {
                update.get(0).tag = cid++;
                center.add(update.get(0));
                continue;
            }
            Node updateNode = update.remove(updateIndex); //Record the nodes that need to be updated
            ArrayList<Integer> maxLabel = new ArrayList<Integer>();//Record the largest number of labels in neighbors
            int[] neighborLabelNum = new int[center.size() + 1];

            // Calculate the largest number of labels among neighbor nodes
            for (int i = 0; i < updateNode.neighbor.size(); i++) {
                if (updateNode.neighbor.get(i).tag != 0) {
                    neighborLabelNum[updateNode.neighbor.get(i).tag]++;
                }
            }
            int maxLabelNum = 0;
            for (int i = 1; i <= center.size(); i++) {
                if (neighborLabelNum[i] > maxLabelNum) {
                    maxLabelNum = neighborLabelNum[i];
                    maxLabel.clear();
                    maxLabel.add(i);
                } else if (neighborLabelNum[i] == maxLabelNum)
                    maxLabel.add(i);
            }
            //If there are multiple candidate tags, select the labelwith the greatest label similarity
            if (maxLabel.size() == 1)
                updateNode.tag = maxLabel.get(0);
            else {
                int simLabelIndex = 0;
                double[] labelsim = new double[center.size() + 1];
                double maxSim = 0;
                for (int i = 0; i < updateNode.neighbor.size(); i++) {
                    if (maxLabel.contains(updateNode.neighbor.get(i).tag)) {
                        labelsim[updateNode.neighbor.get(i).tag] += 1 - jaccard[Integer.parseInt(updateNode.id)][Integer.parseInt(updateNode.neighbor.get(i).id)];
                    }
                }
                for (int i = 1; i <= center.size(); i++) {
                    if (labelsim[i] > maxSim) {
                        maxSim = labelsim[i];
                        simLabelIndex = i;
                    }
                }
                updateNode.tag = simLabelIndex;
            }
        }
        //Free up useless memory space
        jaccard = null;
    }
    public  void readComFile(String comDataSet) throws FileNotFoundException {
        result = new int[node.size()];
        com = new int[node.size() + 1];
        int nodeID = 10000;
        int nodeTag;
        for (int i = 0; i < node.size(); i++)
            result[i] = node.get(i).tag;
        File com_file = new File(comDataSet);
        Scanner in = new Scanner(com_file);
        while (in.hasNext()) {
            nodeID = in.nextInt();
            nodeTag = in.nextInt();
            com[nodeID] = nodeTag;
        }
        if (com[0] == 0) {
            for (int i = 0; i < node.size(); i++)
                com[i] = com[i + 1];
        }
    }

    public  double log2(double x) {
        return (Math.log10(x) / Math.log10(2));
    }

    public  double NMI() {
        int max = 0;
        for (int i = 0; i < node.size(); i++) {
            if (com[i] > max)
                max = com[i];
        }
        c_num = max;
        int[][] m = new int[center.size()][c_num];// 混淆矩阵
        double s1 = 0, s2 = 0, s3 = 0;
        int[] rows = new int[center.size()];
        int[] cols = new int[c_num];
        for (int i = 0; i < node.size(); i++) {
            m[result[i] - 1][com[i] - 1]++;
        }
        for (int i = 0; i < center.size(); i++) {
            for (int j = 0; j < c_num; j++)
                rows[i] = rows[i] + m[i][j];
        }
        for (int i = 0; i < c_num; i++) {
            for (int j = 0; j < center.size(); j++)
                cols[i] = cols[i] + m[j][i];
        }
        for (int i = 0; i < center.size(); i++) {
            for (int j = 0; j < c_num; j++) {
                if (m[i][j] == 0)
                    continue;
                else
                    s1 = s1 + m[i][j] * log2(m[i][j] * node.size() / (double) (rows[i] * cols[j]));
            }
        }

        for (int i = 0; i < center.size(); i++) {
            s2 = s2 + rows[i] * log2((double) rows[i] / node.size());
        }
        for (int i = 0; i < c_num; i++) {
            s3 = s3 + cols[i] * log2((double) cols[i] / node.size());
        }
         return -2*s1/(s2+s3);
    }

    public  double ARI() {
        int[] result_pairs = new int[node.size() * (node.size() - 1) / 2];
        int[] com_pairs = new int[node.size() * (node.size() - 1) / 2];
        int k = 0;
        double a11 = 0, a00 = 0, a01 = 0, a10 = 0;
        for (int i = 0; i < node.size() - 1; i++) {
            for (int j = i + 1; j < node.size(); j++) {
                if (result[i] == result[j])
                    result_pairs[k++] = 1;
                else
                    result_pairs[k++] = 0;
            }
        }
        k = 0;
        for (int i = 0; i < node.size() - 1; i++) {
            for (int j = i + 1; j < node.size(); j++) {
                if (com[i] == com[j])
                    com_pairs[k++] = 1;
                else
                    com_pairs[k++] = 0;
            }
        }
        for (int i = 0; i < node.size() * (node.size() - 1) / 2; i++) {
            if (result_pairs[i] == 1 && com_pairs[i] == 1)
                a11++;
            else if (result_pairs[i] == 0 && com_pairs[i] == 0)
                a00++;
            else if (result_pairs[i] == 0 && com_pairs[i] == 1)
                a10++;
            else if (result_pairs[i] == 1 && com_pairs[i] == 0)
                a01++;
        }

        return (double) 2 * (a00 * a11 - a01 * a10) / ((double) (a00 + a01) * (a01 + a11) + (double)(a00 + a10) * (a10 + a11));
    }

    public  double Q() {
        int[][] ad = new int[node.size()][node.size()];// 邻接矩阵
        int[] du = new int[node.size()];// 存储节点的度
        int[][] membership = new int[node.size()][node.size()];// 节点i,j是否在同一社区，在为1不在为0
        double q = 0;
        result = new int[node.size()];
        for (int i = 0; i < node.size(); i++)
            result[i] = node.get(i).tag;
        for (int i = 0; i < node.size(); i++) {
            for (int j = 0; j < node.size(); j++) {
                if (i == j)
                    membership[i][j] = 1;
                else {
                    if (result[i] == result[j])
                        membership[i][j] = 1;
                    else
                        membership[i][j] = 0;
                }
            }
        }
            for (int i = 0; i < edge.size(); i++) {
                ad[Integer.parseInt(edge.get(i).source)][Integer.parseInt(edge.get(i).target)] = 1;
                ad[Integer.parseInt(edge.get(i).target)][Integer.parseInt(edge.get(i).source)] = 1;
                du[Integer.parseInt(edge.get(i).source)]++;
                du[Integer.parseInt(edge.get(i).target)]++;
            }
            for (int i = 0; i < node.size(); i++) {
                for (int j = 0; j < node.size(); j++) {
                    q = q + (ad[i][j] - du[i] * du[j] / (double) (2 * edge.size())) * membership[i][j];
                }
            }
            q = q / (2 * edge.size());
        return q;
    }
}
