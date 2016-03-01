package org.neo4j.slm;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.stream.Stream;

import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.PerformsWrites;
import org.neo4j.procedure.Procedure;

import static java.lang.String.format;

public class ClusterAllTheThings
{
    @Context
    public org.neo4j.graphdb.GraphDatabaseService db;

    @Procedure
    @PerformsWrites
    public Stream<Cluster> knows(String label) throws IOException
    {
        String query = "MATCH (person1:Person)-[r:KNOWS]->(person2:Person) \n" +
                       "RETURN person1.id AS p1, person2.id AS p2, toFloat(1) AS weight";

        Result rows = db.execute( query );

        ModularityOptimizer.ModularityFunction modularityFunction = ModularityOptimizer.ModularityFunction.Standard;
        Network network = Network.create( modularityFunction, rows );

        double resolution = 1.0;
        int nRandomStarts = 1;
        int nIterations = 10;
        long randomSeed = 0;

        double modularity;

        Random random = new Random( randomSeed );

        double resolution2 = modularityFunction.resolution( resolution, network );

        Map<Integer, Node> cluster = new HashMap<>();
        double maxModularity = Double.NEGATIVE_INFINITY;

        for ( int randomStart = 0; randomStart < nRandomStarts; randomStart++ )
        {
            network.initSingletonClusters();

            int iteration = 0;
            do
            {
                network.runSmartLocalMovingAlgorithm( resolution2, random );
                iteration++;

                modularity = network.calcQualityFunction( resolution2 );
            } while ( (iteration < nIterations) );

            if ( modularity > maxModularity )
            {
                network.orderClustersByNNodes();
                cluster = network.getNodes();
                maxModularity = modularity;
            }
        }

        for ( Map.Entry<Integer, Node> entry : cluster.entrySet() )
        {
            Map<String, Object> params = new HashMap<>();
            params.put("userId", String.valueOf(entry.getKey()));
//            params.put("community", entry.getValue().getCluster());
//            db.execute("MATCH (person:Person {id: {userId}})\n" +
//                       "MERGE (community:Community {id: {communityId}})\n" +
//                       "MERGE (person)-[:IN_COMMUNITY]->(community)",
//                    params);

            try ( Transaction tx = db.beginTx() )
            {
                org.neo4j.graphdb.Node node = db.findNode( Label.label( label ), "id", String.valueOf( entry.getKey() ) );
                node.addLabel( Label.label( (format( "Community-%d`", entry.getValue().getCluster() )) ) );
                tx.success();
            }
        }

        return cluster
                .entrySet()
                .stream()
                .map( ( entry ) -> new Cluster( entry.getKey(), entry.getValue().getCluster() ) );
    }

    public static class Cluster
    {
        public long id;
        public long clusterId;

        public Cluster( int id, int clusterId )
        {
            this.id = id;
            this.clusterId = clusterId;
        }
    }
}
