package org.neo4j.slm;

import org.junit.Test;

import java.io.FileReader;
import java.io.IOException;

import org.neo4j.slm.ModularityOptimizer;
import org.neo4j.slm.Network;

import static org.junit.Assert.*;


public class NetworkTest
{
    @Test
    public void shouldCreateKarateNetwork() throws IOException
    {
        ModularityOptimizer.ModularityFunction alternativeModularityFunction = ModularityOptimizer.ModularityFunction.Standard;
        Network network = Network.create( alternativeModularityFunction, new FileReader( "karate_club_network.txt" ) );

        assertEquals(34, network.getNNodes());

        assertArrayEquals(
                new double[]{16.0, 9.0, 10.0, 6.0, 3.0, 4.0, 4.0, 4.0, 5.0, 2.0, 3.0, 1.0, 2.0, 5.0, 2.0, 2.0, 2.0, 2.0,
                        2.0, 3.0, 2.0, 2.0, 2.0, 5.0, 3.0, 3.0, 2.0, 4.0, 3.0, 4.0, 4.0, 6.0, 12.0, 17.0},
                network.getNodeWeights(), 1.0 );


        assertArrayEquals(
                new double[]{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                        1.0, 1.0, 1.0, 1.0, 1.0},
                network.getEdgeWeights(),
                1.0);
    }
}
