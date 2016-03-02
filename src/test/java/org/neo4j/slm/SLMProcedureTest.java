package org.neo4j.slm;

import org.junit.Test;

import org.neo4j.graphdb.Result;
import org.neo4j.kernel.GraphDatabaseAPI;
import org.neo4j.kernel.impl.proc.Procedures;
import org.neo4j.test.TestGraphDatabaseFactory;

import static org.junit.Assert.assertEquals;

public class SLMProcedureTest
{
    @Test
    public void shouldDetermineClustersGiven0BasedMonotonicIdsWithNoGaps() throws Exception
    {
        GraphDatabaseAPI db = (GraphDatabaseAPI) new TestGraphDatabaseFactory().newImpermanentDatabase();
        db.getDependencyResolver().resolveDependency(Procedures.class).register(SLMProcedure.class);

        // given
        db.execute("CREATE (t1:Topic {id: '0'})-[:SIMILAR {score: 0.33}]->(t2:Topic {id: '1'})").close();

        // when
        Result res = db.execute("CALL org.neo4j.slm.slm('Topic', 'SIMILAR')");

        // then
        assertEquals(0L, res.next().get( "clusterId" ));
        assertEquals(0L, res.next().get( "clusterId" ));
    }

    @Test
    public void shouldDetermineClustersGiven0BasedMonotonicIdsWithGaps() throws Exception
    {
        GraphDatabaseAPI db = (GraphDatabaseAPI) new TestGraphDatabaseFactory().newImpermanentDatabase();
        db.getDependencyResolver().resolveDependency(Procedures.class).register(SLMProcedure.class);

        // given
        db.execute("CREATE (t1:Topic {id: '0'})-[:SIMILAR {score: 0.33}]->(t2:Topic {id: '2'})").close();

        // when
        Result res = db.execute("CALL org.neo4j.slm.slm('Topic', 'SIMILAR')");

        // then
        assertEquals(0L, res.next().get( "clusterId" ));
        assertEquals(0L, res.next().get( "clusterId" ));
    }

}
