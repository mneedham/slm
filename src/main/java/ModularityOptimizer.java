/**
 * ModularityOptimizer
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.2.0, 05/14/14
 */

import java.io.BufferedWriter;
import java.io.Console;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.Random;

public class ModularityOptimizer
{

    enum ModularityFunction
    {
        Standard( 1 )
                {
                    @Override
                    double resolution( double resolution, Network network )
                    {
                        return (resolution / network.getTotalEdgeWeight());
                    }

                    @Override
                    Network createNetwork( int nNodes, int nEdges, int[] firstNeighborIndex, int[] neighbor,
                            double[] edgeWeight, Map<Integer,Node> nodes )
                    {

                        return new Network( nodes );
                    }

                    @Override
                    ReducedNetwork createReducedNetwork( int nNodes, int nEdges, int[] firstNeighborIndex,
                            int[] neighbor,
                            double[] edgeWeight, Map<Integer,Node> nodesMap )
                    {
                        double[] nodeWeight = new double[nNodes];
                        for ( Map.Entry<Integer,Node> entry : nodesMap.entrySet() )
                        {
                            nodeWeight[entry.getKey()] = entry.getValue().weight();
                        }

                        return new ReducedNetwork( nNodes, firstNeighborIndex, neighbor, edgeWeight, nodeWeight, nodesMap );
                    }
                },
        Alternative( 2 )
                {
                    @Override
                    double resolution( double resolution, Network network )
                    {
                        return resolution;
                    }

                    @Override
                    Network createNetwork( int nNodes, int nEdges, int[] firstNeighborIndex, int[] neighbor,
                            double[] edgeWeight2, Map<Integer,Node> nodesMap )
                    {
                        return new Network( neighbor, edgeWeight2 );
                    }

                    @Override
                    ReducedNetwork createReducedNetwork( int nNodes, int nEdges, int[] firstNeighborIndex,
                            int[] neighbor,
                            double[] edgeWeight2, Map<Integer,Node> nodesMap )
                    {
                        return new ReducedNetwork( nNodes, firstNeighborIndex, neighbor, edgeWeight2 );
                    }
                };

        private int id;

        ModularityFunction( int id )
        {
            this.id = id;
        }

        public static ModularityFunction from( int id )
        {
            if ( id == Standard.id )
            {
                return Standard;
            }
            return Alternative;
        }

        abstract double resolution( double resolution, Network network );

        abstract Network createNetwork( int nNodes, int nEdges, int[] firstNeighborIndex, int[] neighbor,
                double[] edgeWeight2, Map<Integer,Node> nodesMap );

        abstract ReducedNetwork createReducedNetwork( int nNodes, int nEdges, int[] firstNeighborIndex, int[] neighbor,
                double[] edgeWeight2, Map<Integer,Node> nodesMap );
    }

    public static void main( String[] args ) throws IOException
    {
        boolean printOutput, update;
        Console console;
        double modularity, maxModularity, resolution;
        int algorithm;
        int i;
        int j;
        ModularityFunction modularityFunction;
        int numberOfClusters;
        int nIterations;
        int nRandomStarts;
        int[] cluster;
        long beginTime, endTime, randomSeed;
        Network network;
        Random random;
        String inputFileName, outputFileName;

        if ( args.length == 9 )
        {
            inputFileName = args[0];
            outputFileName = args[1];
            modularityFunction = ModularityFunction.from( Integer.parseInt( args[2] ) );
            resolution = Double.parseDouble( args[3] );
            algorithm = Integer.parseInt( args[4] );
            nRandomStarts = Integer.parseInt( args[5] );
            nIterations = Integer.parseInt( args[6] );
            randomSeed = Long.parseLong( args[7] );
            printOutput = (Integer.parseInt( args[8] ) > 0);

            if ( printOutput )
            {
                System.out.println( "Modularity Optimizer version 1.2.0 by Ludo Waltman and Nees Jan van Eck" );
                System.out.println();
            }
        }
        else
        {
            console = System.console();
            System.out.println( "Modularity Optimizer version 1.2.0 by Ludo Waltman and Nees Jan van Eck" );
            System.out.println();
            inputFileName = console.readLine( "Input file name: " );
            outputFileName = console.readLine( "Output file name: " );
            modularityFunction =
                    ModularityFunction.from(
                            Integer.parseInt( console.readLine( "Modularity function (1 = standard; 2" +
                                                                " = alternative): " ) ) );
            resolution = Double.parseDouble( console.readLine( "Resolution parameter (e.g., 1.0): " ) );
            algorithm = Integer.parseInt( console.readLine(
                    "Algorithm (1 = Louvain; 2 = Louvain with multilevel refinement; 3 = smart local moving): " ) );
            nRandomStarts = Integer.parseInt( console.readLine( "Number of random starts (e.g., 10): " ) );
            nIterations = Integer.parseInt( console.readLine( "Number of iterations (e.g., 10): " ) );
            randomSeed = Long.parseLong( console.readLine( "Random seed (e.g., 0): " ) );
            printOutput = (Integer.parseInt( console.readLine( "Print output (0 = no; 1 = yes): " ) ) > 0);
            System.out.println();
        }

        if ( printOutput )
        {
            System.out.println( "Reading input file..." );
            System.out.println();
        }

        network = Network.create( inputFileName, modularityFunction );

        if ( printOutput )
        {
            System.out.format( "Number of nodes: %d%n", network.getNNodes() );
            System.out.format( "Number of edges: %d%n", network.getNEdges() / 2 );
            System.out.println();
            System.out.println( "Running " + algorithmDescription( algorithm ) +
                                "..." );
            System.out.println();
        }

        double resolution2 = modularityFunction.resolution( resolution, network );

        beginTime = System.currentTimeMillis();
        cluster = null;
        numberOfClusters = -1;
        maxModularity = Double.NEGATIVE_INFINITY;
        random = new Random( randomSeed );
        for ( i = 0; i < nRandomStarts; i++ )
        {
            if ( printOutput && (nRandomStarts > 1) )
            {
                System.out.format( "Random start: %d%n", i + 1 );
            }

            network.initSingletonClusters();

//            printCurrentClusters( network );

            j = 0;
            update = true;
            do
            {
                if ( printOutput && (nIterations > 1) )
                { System.out.format( "Iteration: %d%n", j + 1 ); }

                if ( algorithm == 1 )
                {
                    update = network.runLouvainAlgorithm( resolution2, random );
                }
                else if ( algorithm == 2 )
                {
                    update = network.runLouvainAlgorithmWithMultilevelRefinement( resolution2, random );
                }
                else if ( algorithm == 3 )
                {
                    network.runSmartLocalMovingAlgorithm( resolution2, random );
//                    printCurrentClusters( network );
                }

                j++;

                modularity = network.calcQualityFunction( resolution2 );

                if ( printOutput && (nIterations > 1) )
                { System.out.format( "Modularity: %.4f%n", modularity ); }
            }
            while ( (j < nIterations) && update );

            if ( modularity > maxModularity )
            {
                network.orderClustersByNNodes();
                cluster = network.getClusters();
                numberOfClusters = network.getNClusters();
                maxModularity = modularity;
            }

            if ( printOutput && (nRandomStarts > 1) )
            {
                if ( nIterations == 1 )
                { System.out.format( "Modularity: %.4f%n", modularity ); }
                System.out.println();
            }
        }
        endTime = System.currentTimeMillis();

        if ( printOutput )
        {
            if ( nRandomStarts == 1 )
            {
                if ( nIterations > 1 )
                { System.out.println(); }
                System.out.format( "Modularity: %.4f%n", maxModularity );
            }
            else
            { System.out.format( "Maximum modularity in %d random starts: %.4f%n", nRandomStarts, maxModularity ); }
            System.out.format( "Number of communities: %d%n", numberOfClusters );
            System.out.format( "Elapsed time: %d seconds%n", Math.round( (endTime - beginTime) / 1000.0 ) );
            System.out.println();
            System.out.println( "Writing output file..." );
            System.out.println();
        }

        writeOutputFile( outputFileName, cluster );
    }

    private static String algorithmDescription( int algorithm )
    {
        return (algorithm == 1) ? "Louvain algorithm" :
               ((algorithm == 2) ? "Louvain algorithm with multilevel refinement" :
                "smart local moving algorithm");
    }

    private static void printCurrentClusters( Network network )
    {
        for ( int item : network.getClusters() )
        {
            System.out.print( item + " " );
        }
        System.out.println();
    }

    private static void writeOutputFile( String fileName, int[] cluster ) throws IOException
    {
        BufferedWriter bufferedWriter;
        int i;

        bufferedWriter = new BufferedWriter( new FileWriter( fileName ) );

        for ( i = 0; i < cluster.length; i++ )
        {
            bufferedWriter.write( Integer.toString( cluster[i] ) );
            bufferedWriter.newLine();
        }

        bufferedWriter.close();
    }
}
