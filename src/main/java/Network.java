/**
 * Network
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.2.0, 05/14/14
 */

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

public class Network implements Cloneable, Serializable
{
    private static final long serialVersionUID = 1;

    private int numberOfNodes;
    private int[] firstNeighborIndex;
    private int[] neighbor;
    private double[] edgeWeight;
    private double totalEdgeWeightSelfLinks;
    private double[] nodeWeight;
    private int nClusters;
    private int[] cluster;

    private double[] clusterWeight;
    private int[] numberNodesPerCluster;
    private int[][] nodePerCluster;
    private boolean clusteringStatsAvailable;

    public int[] getNeighbor()
    {
        return neighbor;
    }

    public static Network load( String fileName ) throws ClassNotFoundException, IOException
    {
        Network network;
        ObjectInputStream objectInputStream;

        objectInputStream = new ObjectInputStream( new FileInputStream( fileName ) );

        network = (Network) objectInputStream.readObject();

        objectInputStream.close();

        return network;
    }

    public Network( int numberOfNodes, int[][] edge )
    {
        this( numberOfNodes, edge, null, null, null );
    }

    public Network( int numberOfNodes, int[][] edge, double[] edgeWeight )
    {
        this( numberOfNodes, edge, edgeWeight, null, null );
    }

    public Network( int numberOfNodes, int[][] edge, double[] edgeWeight, double[] nodeWeight )
    {
        this( numberOfNodes, edge, edgeWeight, nodeWeight, null );
    }

    public Network( int numberOfNodes, int[][] edge, double[] edgeWeight, double[] nodeWeight, int[] cluster )
    {
        double[] edgeWeight2;
        int i, j, nEdges, nEdgesWithoutSelfLinks;
        int[] neighbor;

        this.numberOfNodes = numberOfNodes;

        nEdges = edge[0].length;
        firstNeighborIndex = new int[numberOfNodes + 1];
        if ( edgeWeight == null )
        {
            edgeWeight = new double[nEdges];
            for ( i = 0; i < nEdges; i++ )
            { edgeWeight[i] = 1; }
        }
        totalEdgeWeightSelfLinks = 0;

        neighbor = new int[nEdges];
        edgeWeight2 = new double[nEdges];
        i = 1;
        nEdgesWithoutSelfLinks = 0;
        for ( j = 0; j < nEdges; j++ )
        {
            if ( edge[0][j] == edge[1][j] )
            { totalEdgeWeightSelfLinks += edgeWeight[j]; }
            else
            {
                if ( edge[0][j] >= i )
                {
                    for (; i <= edge[0][j]; i++ )
                    { firstNeighborIndex[i] = nEdgesWithoutSelfLinks; }
                }
                neighbor[nEdgesWithoutSelfLinks] = edge[1][j];
                edgeWeight2[nEdgesWithoutSelfLinks] = edgeWeight[j];
                nEdgesWithoutSelfLinks++;
            }
        }
        for (; i <= numberOfNodes; i++ )
        { firstNeighborIndex[i] = nEdgesWithoutSelfLinks; }

        this.neighbor = new int[nEdgesWithoutSelfLinks];
        System.arraycopy( neighbor, 0, this.neighbor, 0, nEdgesWithoutSelfLinks );
        this.edgeWeight = new double[nEdgesWithoutSelfLinks];
        System.arraycopy( edgeWeight2, 0, this.edgeWeight, 0, nEdgesWithoutSelfLinks );

        if ( nodeWeight == null )
        {
            this.nodeWeight = new double[numberOfNodes];
            for ( i = 0; i < numberOfNodes; i++ )
            { this.nodeWeight[i] = 1; }
        }
        else
        { this.nodeWeight = nodeWeight; }

        setClusters( cluster );
    }

    public Network( int numberOfNodes, int[] firstNeighborIndex, int[] neighbor )
    {
        this( numberOfNodes, firstNeighborIndex, neighbor, null, null, null );
    }

    public Network( int numberOfNodes, int[] firstNeighborIndex, int[] neighbor, double[] edgeWeight )
    {
        this( numberOfNodes, firstNeighborIndex, neighbor, edgeWeight, null, null );
    }

    public Network( int numberOfNodes, int[] firstNeighborIndex, int[] neighbor, double[] edgeWeight, double[] nodeWeight )
    {
        this( numberOfNodes, firstNeighborIndex, neighbor, edgeWeight, nodeWeight, null );
    }

    public Network( int numberOfNodes, int[] firstNeighborIndex, int[] neighbor, double[] edgeWeight, double[] nodeWeight,
            int[] cluster )
    {
        int i, nEdges;

        this.numberOfNodes = numberOfNodes;

        this.firstNeighborIndex = firstNeighborIndex;
        this.neighbor = neighbor;

        if ( edgeWeight == null )
        {
            nEdges = neighbor.length;
            this.edgeWeight = new double[nEdges];
            for ( i = 0; i < nEdges; i++ )
            { this.edgeWeight[i] = 1; }
        }
        else
        { this.edgeWeight = edgeWeight; }

        if ( nodeWeight == null )
        {
            this.nodeWeight = new double[numberOfNodes];
            for ( i = 0; i < numberOfNodes; i++ )
            { this.nodeWeight[i] = 1; }
        }
        else
        { this.nodeWeight = nodeWeight; }

        setClusters( cluster );
    }

    public static Network create( String fileName, ModularityOptimizer.ModularityFunction modularityFunction ) throws IOException
    {
        double[] edgeWeight1, edgeWeight2;
        int i, j, nEdges, nNodes;
        int[] neighbor, source, destination;

        int numberOfLines = numberOfLines( fileName );
        try ( BufferedReader bufferedReader = new BufferedReader( new FileReader( fileName ) ) )
        {
            source = new int[numberOfLines];
            destination = new int[numberOfLines];
            edgeWeight1 = new double[numberOfLines];

            Set<Integer> nodes = new HashSet<>(  );
            for ( j = 0; j < numberOfLines; j++ )
            {
                String[] splittedLine = bufferedReader.readLine().split( "\t" );

                source[j] = Integer.parseInt( splittedLine[0] );
                destination[j] = Integer.parseInt( splittedLine[1] );

                nodes.add( source[j] );
                nodes.add( destination[j] );

                edgeWeight1[j] = (splittedLine.length > 2) ? Double.parseDouble( splittedLine[2] ) : 1;
            }
            nNodes = nodes.size();
        }

        int[] numberOfNeighbours = numberOfNeighbours( nNodes, source, destination, numberOfLines );
        int[] firstNeighborIndex = new int[nNodes + 1];
        nEdges = 0;
        for ( i = 0; i < nNodes; i++ )
        {
            firstNeighborIndex[i] = nEdges;
            nEdges += numberOfNeighbours[i];
        }
        firstNeighborIndex[nNodes] = nEdges;

        neighbor = new int[nEdges];
        edgeWeight2 = new double[nEdges];

        Arrays.fill( numberOfNeighbours, 0 );
        for ( i = 0; i < numberOfLines; i++ )
        {
            if ( source[i] < destination[i] )
            {
                j = firstNeighborIndex[source[i]] + numberOfNeighbours[source[i]];
                neighbor[j] = destination[i];
                edgeWeight2[j] = edgeWeight1[i];
                numberOfNeighbours[source[i]]++;

                j = firstNeighborIndex[destination[i]] + numberOfNeighbours[destination[i]];
                neighbor[j] = source[i];
                edgeWeight2[j] = edgeWeight1[i];
                numberOfNeighbours[destination[i]]++;
            }
        }

        return modularityFunction.createNetwork(nNodes, nEdges, firstNeighborIndex, neighbor, edgeWeight2 );
    }

    private static int[] numberOfNeighbours( int nNodes, int[] source, int[] destination, int numberOfLines )
    {
        int i;
        int[] numberOfNeighbours = new int[nNodes];
        for ( i = 0; i < numberOfLines; i++ )
        {
            if ( source[i] < destination[i] )
            {
                numberOfNeighbours[source[i]]++;
                numberOfNeighbours[destination[i]]++;
            }
        }
        return numberOfNeighbours;
    }

    private static int numberOfLines( String fileName ) throws IOException
    {
        BufferedReader bufferedReader = new BufferedReader( new FileReader( fileName ) );
        int nLines;

        nLines = 0;
        while ( bufferedReader.readLine() != null )
        { nLines++; }

        bufferedReader.close();
        return nLines;
    }

    public Object clone()
    {
        Network clonedNetwork;

        try
        {
            clonedNetwork = (Network) super.clone();

            if ( cluster != null )
            { clonedNetwork.cluster = (int[]) cluster.clone(); }
            clonedNetwork.deleteClusteringStats();

            return clonedNetwork;
        }
        catch ( CloneNotSupportedException e )
        {
            return null;
        }
    }

    public void save( String fileName ) throws IOException
    {
        ObjectOutputStream objectOutputStream;

        objectOutputStream = new ObjectOutputStream( new FileOutputStream( fileName ) );

        objectOutputStream.writeObject( this );

        objectOutputStream.close();
    }

    public int getNNodes()
    {
        return numberOfNodes;
    }

    public int getNEdges()
    {
        return neighbor.length;
    }

    public int[][] getEdges()
    {
        int[][] edge;
        int i, j;

        edge = new int[2][neighbor.length];
        for ( i = 0; i < numberOfNodes; i++ )
        {
            for ( j = firstNeighborIndex[i]; j < firstNeighborIndex[i + 1]; j++ )
            {
                edge[0][j] = i;
                edge[1][j] = neighbor[j];
            }
        }

        return edge;
    }

    public double getTotalEdgeWeight()
    {
        double totalEdgeWeight;
        int i;

        totalEdgeWeight = totalEdgeWeightSelfLinks;
        for ( i = 0; i < neighbor.length; i++ )
        {
            totalEdgeWeight += edgeWeight[i];
        }

        return totalEdgeWeight;
    }

    public double[] getEdgeWeights()
    {
        return edgeWeight;
    }

    public double getTotalNodeWeight()
    {
        double totalNodeWeight;
        int i;

        totalNodeWeight = 0;
        for ( i = 0; i < numberOfNodes; i++ )
        { totalNodeWeight += nodeWeight[i]; }

        return totalNodeWeight;
    }

    public double[] getNodeWeights()
    {
        return nodeWeight;
    }

    public int getNClusters()
    {
        return nClusters;
    }

    public int[] getClusters()
    {
        return cluster;
    }

    public double[] getClusterWeights()
    {
        if ( cluster == null )
        { return null; }

        if ( !clusteringStatsAvailable )
        { calcClusteringStats(); }

        return clusterWeight;
    }

    public int[] getNNodesPerCluster()
    {
        if ( cluster == null )
        { return null; }

        if ( !clusteringStatsAvailable )
        { calcClusteringStats(); }

        return numberNodesPerCluster;
    }

    public int[][] getNodesPerCluster()
    {
        if ( cluster == null )
        { return null; }

        if ( !clusteringStatsAvailable )
        { calcClusteringStats(); }

        return nodePerCluster;
    }

    public void setClusters( int[] cluster )
    {
        nClusters = Clusters.calculateNumberOfClusters( cluster, numberOfNodes );
        this.cluster = cluster;
        deleteClusteringStats();
    }

    public void initSingletonClusters()
    {
        int i;

        nClusters = numberOfNodes;
        cluster = new int[numberOfNodes];
        for ( i = 0; i < numberOfNodes; i++ )
        {
            cluster[i] = i;
        }

        deleteClusteringStats();
    }

    public void findConnectedComponents()
    {
        int i, j;
        int[] neighborIndex, node;

        cluster = new int[numberOfNodes];
        init();

        node = new int[numberOfNodes];
        neighborIndex = new int[numberOfNodes];

        nClusters = 0;
        for ( i = 0; i < numberOfNodes; i++ )
        {
            if ( cluster[i] == -1 )
            {
                cluster[i] = nClusters;
                node[0] = i;
                neighborIndex[0] = firstNeighborIndex[i];
                j = 0;
                do
                {
                    if ( neighborIndex[j] == firstNeighborIndex[node[j] + 1] )
                    { j--; }
                    else if ( cluster[neighbor[neighborIndex[j]]] == -1 )
                    {
                        cluster[neighbor[neighborIndex[j]]] = nClusters;
                        node[j + 1] = neighbor[neighborIndex[j]];
                        neighborIndex[j + 1] = firstNeighborIndex[node[j + 1]];
                        neighborIndex[j]++;
                        j++;
                    }
                    else
                    { neighborIndex[j]++; }
                }
                while ( j >= 0 );

                nClusters++;
            }
        }

        deleteClusteringStats();
    }

    private void init()
    {
        int i;
        for ( i = 0; i < numberOfNodes; i++ )
        { cluster[i] = -1; }
    }

    private static void printCluster( int[] clusters )
    {
        for ( int item : clusters )
        {
            System.out.print( item + " " );
        }
        System.out.println();
    }

    public void mergeClusters( int[] newCluster )
    {
        if ( cluster == null )
        {
            return;
        }

        int i = 0;
        for ( int j = 0; j < numberOfNodes; j++ )
        {
            int k = newCluster[cluster[j]];
            if ( k > i )
            {
                i = k;
            }
            cluster[j] = k;
        }
        nClusters = i + 1;

        deleteClusteringStats();
    }

    public boolean removeCluster( int cluster )
    {
        boolean removed;

        if ( this.cluster == null )
        { return false; }

        if ( !clusteringStatsAvailable )
        { calcClusteringStats(); }

        removed = removeCluster2( cluster );

        deleteClusteringStats();

        return removed;
    }

    public void removeSmallClusters( double minClusterWeight )
    {
        boolean[] ignore;
        double minClusterWeight2;
        int i, smallestCluster;
        Network reducedNetwork;

        if ( cluster == null )
        { return; }

        reducedNetwork = getReducedNetwork();
        reducedNetwork.initSingletonClusters();
        reducedNetwork.calcClusteringStats();

        ignore = new boolean[nClusters];
        do
        {
            smallestCluster = -1;
            minClusterWeight2 = minClusterWeight;
            for ( i = 0; i < reducedNetwork.nClusters; i++ )
            {
                if ( (!ignore[i]) && (reducedNetwork.clusterWeight[i] < minClusterWeight2) )
                {
                    smallestCluster = i;
                    minClusterWeight2 = reducedNetwork.clusterWeight[i];
                }
            }

            if ( smallestCluster >= 0 )
            {
                reducedNetwork.removeCluster2( smallestCluster );
                ignore[smallestCluster] = true;
            }
        }
        while ( smallestCluster >= 0 );

        mergeClusters( reducedNetwork.getClusters() );
    }

    public void orderClustersByWeight()
    {
        orderClusters( true );
    }

    public void orderClustersByNNodes()
    {
        orderClusters( false );
    }

    public Network getSubnetwork( int cluster )
    {
        double[] subnetworkEdgeWeight;
        int[] subnetworkNeighbor, subnetworkNode;

        if ( this.cluster == null )
        { return null; }

        if ( !clusteringStatsAvailable )
        { calcClusteringStats(); }

        subnetworkNode = new int[numberOfNodes];
        subnetworkNeighbor = new int[neighbor.length];
        subnetworkEdgeWeight = new double[edgeWeight.length];

        return getSubnetwork( cluster, subnetworkNode, subnetworkNeighbor, subnetworkEdgeWeight );
    }

    public Network[] getSubnetworks()
    {
        if ( cluster == null )
        {
            return null;
        }

        if ( !clusteringStatsAvailable )
        {
            calcClusteringStats();
        }

        Network[] subnetwork = new Network[nClusters];
        int[] subnetworkNode = new int[numberOfNodes];
        int[] subnetworkNeighbor = new int[neighbor.length];
        double[] subnetworkEdgeWeight = new double[edgeWeight.length];

        for ( int i = 0; i < nClusters; i++ )
        {
            subnetwork[i] = getSubnetwork( i, subnetworkNode, subnetworkNeighbor, subnetworkEdgeWeight );
        }

        return subnetwork;
    }

    public Network getReducedNetwork()
    {
        double[] reducedNetworkEdgeWeight1, reducedNetworkEdgeWeight2;
        int i, j, k, l, m, reducedNetworkNEdges1, reducedNetworkNEdges2;
        int[] reducedNetworkNeighbor1, reducedNetworkNeighbor2;
        Network reducedNetwork;

        if ( cluster == null )
        { return null; }

        if ( !clusteringStatsAvailable )
        { calcClusteringStats(); }

        reducedNetwork = new Network();

        reducedNetwork.numberOfNodes = nClusters;
        reducedNetwork.firstNeighborIndex = new int[nClusters + 1];
        reducedNetwork.totalEdgeWeightSelfLinks = totalEdgeWeightSelfLinks;
        reducedNetwork.nodeWeight = new double[nClusters];

        reducedNetworkNeighbor1 = new int[neighbor.length];
        reducedNetworkEdgeWeight1 = new double[edgeWeight.length];

        reducedNetworkNeighbor2 = new int[nClusters - 1];
        reducedNetworkEdgeWeight2 = new double[nClusters];

        reducedNetworkNEdges1 = 0;
        for ( i = 0; i < nClusters; i++ )
        {
            reducedNetworkNEdges2 = 0;
            for ( j = 0; j < nodePerCluster[i].length; j++ )
            {
                k = nodePerCluster[i][j];

                for ( l = firstNeighborIndex[k]; l < firstNeighborIndex[k + 1]; l++ )
                {
                    m = cluster[neighbor[l]];
                    if ( m != i )
                    {
                        if ( reducedNetworkEdgeWeight2[m] == 0 )
                        {
                            reducedNetworkNeighbor2[reducedNetworkNEdges2] = m;
                            reducedNetworkNEdges2++;
                        }
                        reducedNetworkEdgeWeight2[m] += edgeWeight[l];
                    }
                    else
                    {
                        reducedNetwork.totalEdgeWeightSelfLinks += edgeWeight[l];
                    }
                }

                reducedNetwork.nodeWeight[i] += nodeWeight[k];
            }

            for ( j = 0; j < reducedNetworkNEdges2; j++ )
            {
                reducedNetworkNeighbor1[reducedNetworkNEdges1 + j] = reducedNetworkNeighbor2[j];
                reducedNetworkEdgeWeight1[reducedNetworkNEdges1 + j] = reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[j]];
                reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[j]] = 0;
            }
            reducedNetworkNEdges1 += reducedNetworkNEdges2;

            reducedNetwork.firstNeighborIndex[i + 1] = reducedNetworkNEdges1;
        }

        reducedNetwork.neighbor = new int[reducedNetworkNEdges1];
        reducedNetwork.edgeWeight = new double[reducedNetworkNEdges1];
        System.arraycopy( reducedNetworkNeighbor1, 0, reducedNetwork.neighbor, 0, reducedNetworkNEdges1 );
        System.arraycopy( reducedNetworkEdgeWeight1, 0, reducedNetwork.edgeWeight, 0, reducedNetworkNEdges1 );

        return reducedNetwork;
    }

    public Network getLargestConnectedComponent()
    {
        double maxClusterWeight;
        int i, largestCluster;
        Network clonedNetwork;

        clonedNetwork = (Network) clone();

        clonedNetwork.findConnectedComponents();

        clonedNetwork.calcClusteringStats();
        largestCluster = -1;
        maxClusterWeight = -1;
        for ( i = 0; i < clonedNetwork.nClusters; i++ )
        {
            if ( clonedNetwork.clusterWeight[i] > maxClusterWeight )
            {
                largestCluster = i;
                maxClusterWeight = clonedNetwork.clusterWeight[i];
            }
        }

        return clonedNetwork.getSubnetwork( largestCluster );
    }

    public double calcQualityFunction( double resolution )
    {
        double qualityFunction, totalEdgeWeight;
        int i, j, k;

        if ( cluster == null )
        { return Double.NaN; }

        if ( !clusteringStatsAvailable )
        { calcClusteringStats(); }

        qualityFunction = totalEdgeWeightSelfLinks;
        totalEdgeWeight = totalEdgeWeightSelfLinks;
        for ( i = 0; i < numberOfNodes; i++ )
        {
            j = cluster[i];
            for ( k = firstNeighborIndex[i]; k < firstNeighborIndex[i + 1]; k++ )
            {
                if ( cluster[neighbor[k]] == j )
                { qualityFunction += edgeWeight[k]; }
                totalEdgeWeight += edgeWeight[k];
            }
        }

        for ( i = 0; i < nClusters; i++ )
        { qualityFunction -= clusterWeight[i] * clusterWeight[i] * resolution; }

        qualityFunction /= totalEdgeWeight;

        return qualityFunction;
    }

    public boolean runLocalMovingAlgorithm( double resolution )
    {
        return runLocalMovingAlgorithm( resolution, new Random() );
    }

    public boolean runLocalMovingAlgorithm( double resolution, Random random )
    {
        double qualityFunction;
        int[] newCluster;

        if ( (cluster == null) || (numberOfNodes == 1) )
        {
            return false;
        }

        boolean update = false;
        double[] clusterWeight = calculateClusterWeight( numberOfNodes, cluster, nodeWeight );
        int[] numberOfNodesPerCluster = calculateNumberOfNodesPerCluster( numberOfNodes, cluster );

        int numberUnusedClusters = 0;
        int[] unusedCluster = new int[numberOfNodes];
        for ( int i = 0; i < numberOfNodes; i++ )
        {
            if ( numberOfNodesPerCluster[i] == 0 )
            {
                unusedCluster[numberUnusedClusters] = i;
                numberUnusedClusters++;
            }
        }

        int[] nodeOrder = nodeOrder( numberOfNodes, random );
        double[] edgeWeightPerCluster = new double[numberOfNodes];
        int[] neighboringCluster = new int[numberOfNodes - 1];

        int numberStableNodes = 0;
        int i = 0;
        do
        {
            int j = nodeOrder[i];

            int numberOfNeighbouringClusters = 0;
            for ( int k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++ )
            {
                int l = cluster[neighbor[k]];
                if ( edgeWeightPerCluster[l] == 0 )
                {
                    neighboringCluster[numberOfNeighbouringClusters] = l;
                    numberOfNeighbouringClusters++;
                }
                edgeWeightPerCluster[l] += edgeWeight[k];
            }

            clusterWeight[cluster[j]] -= nodeWeight[j];
            numberOfNodesPerCluster[cluster[j]]--;
            if ( numberOfNodesPerCluster[cluster[j]] == 0 )
            {
                unusedCluster[numberUnusedClusters] = cluster[j];
                numberUnusedClusters++;
            }

            int bestCluster = -1;
            double maxQualityFunction = 0;
            for ( int k = 0; k < numberOfNeighbouringClusters; k++ )
            {
                int l = neighboringCluster[k];
                qualityFunction = edgeWeightPerCluster[l] - nodeWeight[j] * clusterWeight[l] * resolution;
                if ( (qualityFunction > maxQualityFunction) ||
                     ((qualityFunction == maxQualityFunction) && (l < bestCluster)) )
                {
                    bestCluster = l;
                    maxQualityFunction = qualityFunction;
                }
                edgeWeightPerCluster[l] = 0;
            }
            if ( maxQualityFunction == 0 )
            {
                bestCluster = unusedCluster[numberUnusedClusters - 1];
                numberUnusedClusters--;
            }

            clusterWeight[bestCluster] += nodeWeight[j];
            numberOfNodesPerCluster[bestCluster]++;
            if ( bestCluster == cluster[j] )
            {
                numberStableNodes++;
            }
            else
            {
                cluster[j] = bestCluster;
                numberStableNodes = 1;
                update = true;
            }

            i = (i < numberOfNodes - 1) ? (i + 1) : 0;
        } while ( numberStableNodes < numberOfNodes );

        newCluster = new int[numberOfNodes];
        nClusters = 0;
        for ( i = 0; i < numberOfNodes; i++ )
        {
            if ( numberOfNodesPerCluster[i] > 0 )
            {
                newCluster[i] = nClusters;
                nClusters++;
            }
        }
        for ( i = 0; i < numberOfNodes; i++ )
        {
            cluster[i] = newCluster[cluster[i]];
        }

        deleteClusteringStats();

        return update;
    }

    private int[] nodeOrder( int numberOfNodes, Random random )
    {
        int[] nodeOrder = new int[numberOfNodes];
        for ( int i = 0; i < numberOfNodes; i++ )
        {
            nodeOrder[i] = i;
        }

        for ( int i = 0; i < numberOfNodes; i++ )
        {
            int j = random.nextInt( numberOfNodes );
            int k = nodeOrder[i];
            nodeOrder[i] = nodeOrder[j];
            nodeOrder[j] = k;
        }
        return nodeOrder;
    }

    private int[] calculateNumberOfNodesPerCluster( int numberOfNodes, int[] cluster )
    {
        int[] numberOfNodesPerCluster = new int[numberOfNodes];
        for ( int i = 0; i < numberOfNodes; i++ )
        {
            numberOfNodesPerCluster[cluster[i]]++;
        }
        return numberOfNodesPerCluster;
    }

    private double[] calculateClusterWeight( int nNodes, int[] cluster, double[] nodeWeight )
    {
        double[] clusterWeight = new double[nNodes];
        for ( int i = 0; i < nNodes; i++ )
        {
            clusterWeight[cluster[i]] += nodeWeight[i];
        }
        return clusterWeight;
    }

    public boolean runLouvainAlgorithm( double resolution )
    {
        return runLouvainAlgorithm( resolution, new Random() );
    }

    public boolean runLouvainAlgorithm( double resolution, Random random )
    {
        boolean update, update2;
        Network reducedNetwork;

        if ( (cluster == null) || (numberOfNodes == 1) )
        { return false; }

        update = runLocalMovingAlgorithm( resolution, random );

        if ( nClusters < numberOfNodes )
        {
            reducedNetwork = getReducedNetwork();
            reducedNetwork.initSingletonClusters();

            update2 = reducedNetwork.runLouvainAlgorithm( resolution, random );

            if ( update2 )
            {
                update = true;

                mergeClusters( reducedNetwork.getClusters() );
            }
        }

        deleteClusteringStats();

        return update;
    }

    public boolean runLouvainAlgorithmWithMultilevelRefinement( double resolution )
    {
        return runLouvainAlgorithmWithMultilevelRefinement( resolution, new Random() );
    }

    public boolean runLouvainAlgorithmWithMultilevelRefinement( double resolution, Random random )
    {
        boolean update, update2;
        Network reducedNetwork;

        if ( (cluster == null) || (numberOfNodes == 1) )
        { return false; }

        update = runLocalMovingAlgorithm( resolution, random );

        if ( nClusters < numberOfNodes )
        {
            reducedNetwork = getReducedNetwork();
            reducedNetwork.initSingletonClusters();

            update2 = reducedNetwork.runLouvainAlgorithm( resolution, random );

            if ( update2 )
            {
                update = true;

                mergeClusters( reducedNetwork.getClusters() );

                runLocalMovingAlgorithm( resolution, random );
            }
        }

        deleteClusteringStats();

        return update;
    }

    public boolean runSmartLocalMovingAlgorithm( double resolution )
    {
        return runSmartLocalMovingAlgorithm( resolution, new Random() );
    }

    public boolean runSmartLocalMovingAlgorithm( double resolution, Random random )
    {
        int i, j, k;
        int[] reducedNetworkCluster, subnetworkCluster;
        Network reducedNetwork;

        if ( (cluster == null) || (numberOfNodes == 1) )
        {
            return false;
        }

        boolean update = runLocalMovingAlgorithm( resolution, random );
        // after this cluster is updated with the 8 clusters

        if ( nClusters < numberOfNodes )
        {
            if ( !clusteringStatsAvailable )
            {
                calcClusteringStats();
            }

            Network[] subnetwork = getSubnetworks();

            nClusters = 0;
            for ( i = 0; i < subnetwork.length; i++ )
            {
                subnetwork[i].initSingletonClusters();
                subnetwork[i].runLocalMovingAlgorithm( resolution, random );

                subnetworkCluster = subnetwork[i].getClusters();
                for ( j = 0; j < subnetworkCluster.length; j++ )
                {
                    cluster[nodePerCluster[i][j]] = nClusters + subnetworkCluster[j];
                }
                nClusters += subnetwork[i].getNClusters();
            }
            calcClusteringStats();

            reducedNetwork = getReducedNetwork();

            reducedNetworkCluster = new int[nClusters];
            i = 0;
            for ( j = 0; j < subnetwork.length; j++ )
            {
                for ( k = 0; k < subnetwork[j].getNClusters(); k++ )
                {
                    reducedNetworkCluster[i] = j;
                    i++;
                }
            }
            reducedNetwork.setClusters( reducedNetworkCluster );

            update |= reducedNetwork.runSmartLocalMovingAlgorithm( resolution, random );

            mergeClusters( reducedNetwork.getClusters() );
//            printCluster( cluster );
        }

        deleteClusteringStats();

        return update;
    }

    private Network()
    {
    }

    private void writeObject( ObjectOutputStream out ) throws IOException
    {
        deleteClusteringStats();

        out.defaultWriteObject();
    }

    private boolean removeCluster2( int cluster )
    {
        double maxQualityFunction, qualityFunction;
        double[] reducedNetworkEdgeWeight;
        int bestCluster, i, j;

        reducedNetworkEdgeWeight = new double[nClusters];
        for ( i = 0; i < numberOfNodes; i++ )
        {
            if ( this.cluster[i] == cluster )
            {
                for ( j = firstNeighborIndex[i]; j < firstNeighborIndex[i + 1]; j++ )
                { reducedNetworkEdgeWeight[this.cluster[neighbor[j]]] += edgeWeight[j]; }
            }
        }

        bestCluster = -1;
        maxQualityFunction = 0;
        for ( i = 0; i < nClusters; i++ )
        {
            if ( (i != cluster) && (clusterWeight[i] > 0) )
            {
                qualityFunction = reducedNetworkEdgeWeight[i] / clusterWeight[i];
                if ( qualityFunction > maxQualityFunction )
                {
                    bestCluster = i;
                    maxQualityFunction = qualityFunction;
                }
            }
        }

        if ( bestCluster == -1 )
        { return false; }

        for ( i = 0; i < numberOfNodes; i++ )
        {
            if ( this.cluster[i] == cluster )
            { this.cluster[i] = bestCluster; }
        }

        clusterWeight[bestCluster] += clusterWeight[cluster];
        clusterWeight[cluster] = 0;

        if ( cluster == nClusters - 1 )
        {
            i = 0;
            for ( j = 0; j < numberOfNodes; j++ )
            {
                if ( this.cluster[j] > i )
                { i = this.cluster[j]; }
            }
            nClusters = i + 1;
        }

        return true;
    }

    private void orderClusters( boolean orderByWeight )
    {
        class ClusterSize implements Comparable<ClusterSize>
        {
            public int cluster;
            public double size;

            public ClusterSize( int cluster, double size )
            {
                this.cluster = cluster;
                this.size = size;
            }

            public int compareTo( ClusterSize cluster )
            {
                return (cluster.size > size) ? 1 : ((cluster.size < size) ? -1 : 0);
            }
        }

        ClusterSize[] clusterSize;
        int i;
        int[] newCluster;

        if ( cluster == null )
        { return; }

        if ( !clusteringStatsAvailable )
        { calcClusteringStats(); }

        clusterSize = new ClusterSize[nClusters];
        for ( i = 0; i < nClusters; i++ )
        { clusterSize[i] = new ClusterSize( i, orderByWeight ? clusterWeight[i] : numberNodesPerCluster[i] ); }

        Arrays.sort( clusterSize );

        newCluster = new int[nClusters];
        i = 0;
        do
        {
            newCluster[clusterSize[i].cluster] = i;
            i++;
        }
        while ( (i < nClusters) && (clusterSize[i].size > 0) );
        nClusters = i;
        for ( i = 0; i < numberOfNodes; i++ )
        { cluster[i] = newCluster[cluster[i]]; }

        deleteClusteringStats();
    }

    private Network getSubnetwork( int cluster, int[] subnetworkNode, int[] subnetworkNeighbor,
            double[] subnetworkEdgeWeight )
    {
        int i, j, k;

        Network subnetwork = new Network();

        int subnetworkNNodes = nodePerCluster[cluster].length;
        subnetwork.numberOfNodes = subnetworkNNodes;

        if ( subnetworkNNodes == 1 )
        {
            subnetwork.firstNeighborIndex = new int[2];
            subnetwork.neighbor = new int[0];
            subnetwork.edgeWeight = new double[0];
            subnetwork.nodeWeight = new double[]{nodeWeight[nodePerCluster[cluster][0]]};
        }
        else
        {
            for ( i = 0; i < nodePerCluster[cluster].length; i++ )
            {
                subnetworkNode[nodePerCluster[cluster][i]] = i;
            }

            subnetwork.firstNeighborIndex = new int[subnetworkNNodes + 1];
            subnetwork.nodeWeight = new double[subnetworkNNodes];

            int subnetworkNEdges = 0;
            for ( i = 0; i < subnetworkNNodes; i++ )
            {
                j = nodePerCluster[cluster][i];

                for ( k = firstNeighborIndex[j]; k < firstNeighborIndex[j + 1]; k++ )
                {
                    if ( this.cluster[neighbor[k]] == cluster )
                    {
                        subnetworkNeighbor[subnetworkNEdges] = subnetworkNode[neighbor[k]];
                        subnetworkEdgeWeight[subnetworkNEdges] = edgeWeight[k];
                        subnetworkNEdges++;
                    }
                }

                subnetwork.firstNeighborIndex[i + 1] = subnetworkNEdges;
                subnetwork.nodeWeight[i] = nodeWeight[j];
            }

            subnetwork.neighbor = new int[subnetworkNEdges];
            subnetwork.edgeWeight = new double[subnetworkNEdges];
            System.arraycopy( subnetworkNeighbor, 0, subnetwork.neighbor, 0, subnetworkNEdges );
            System.arraycopy( subnetworkEdgeWeight, 0, subnetwork.edgeWeight, 0, subnetworkNEdges );
        }

        subnetwork.totalEdgeWeightSelfLinks = 0;

        return subnetwork;
    }

    private void calcClusteringStats()
    {
        int i, j;

        clusterWeight = new double[nClusters];
        numberNodesPerCluster = new int[nClusters];
        nodePerCluster = new int[nClusters][];

        for ( i = 0; i < numberOfNodes; i++ )
        {
            clusterWeight[cluster[i]] += nodeWeight[i];
            numberNodesPerCluster[cluster[i]]++;
        }

        for ( i = 0; i < nClusters; i++ )
        {
            nodePerCluster[i] = new int[numberNodesPerCluster[i]];
            numberNodesPerCluster[i] = 0;
        }

        for ( i = 0; i < numberOfNodes; i++ )
        {
            j = cluster[i];
            nodePerCluster[j][numberNodesPerCluster[j]] = i;
            numberNodesPerCluster[j]++;
        }

        clusteringStatsAvailable = true;
    }

    private void deleteClusteringStats()
    {
        clusterWeight = null;
        numberNodesPerCluster = null;
        nodePerCluster = null;

        clusteringStatsAvailable = false;
    }
}
