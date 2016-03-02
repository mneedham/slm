package org.neo4j.slm;

//import java.util.Map;
//
//import org.junit.Test;
//
//import org.neo4j.graphdb.GraphDatabaseService;
//import org.neo4j.graphdb.Result;
//import org.neo4j.graphdb.factory.GraphDatabaseFactory;
//import org.neo4j.kernel.GraphDatabaseAPI;
//import org.neo4j.kernel.impl.proc.Procedures;
//
//import static junit.framework.TestCase.assertFalse;
//import static org.junit.Assert.assertEquals;
//import static org.junit.Assert.assertTrue;

import java.util.Map;

import org.junit.Test;

import org.neo4j.graphdb.Result;
import org.neo4j.kernel.GraphDatabaseAPI;
import org.neo4j.kernel.impl.proc.Procedures;
import org.neo4j.test.TestGraphDatabaseFactory;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class SLMProcedureTest
{
    @Test
    public void shouldBeEpic() throws Exception
    {
        GraphDatabaseAPI db = (GraphDatabaseAPI) new TestGraphDatabaseFactory().newImpermanentDatabase();

        db.getDependencyResolver().resolveDependency(Procedures.class).register(SLMProcedure.class);

        // given Alice knowing Bob and Charlie and Dan knowing no-one
        db.execute("CREATE (t1:Topic {id: '1'})-[:SIMILAR {score: 0.33}]->(t2:Topic {id: '2'})").close();

        // when retrieving the degree of the User label
        Result res = db.execute("CALL org.neo4j.slm.slm('Topic', 'SIMILAR')");

        // then we expect one result-row with min-degree 0 and max-degree 2
        assertTrue(res.hasNext());
        Map<String,Object> row = res.next();
        assertEquals("User", row.get("label"));
        assertEquals(0L, row.get("min"));
        assertEquals(2L, row.get("max"));
        assertEquals(4L, row.get("count"));
        assertFalse(res.hasNext());
    }


}
