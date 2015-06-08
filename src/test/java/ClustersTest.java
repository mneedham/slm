import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by markneedham on 07/06/15.
 */
public class ClustersTest
{
    @Test
    public void shouldCalculateNumberOfClusters()
    {

        assertEquals( 9,
                Clusters.calculateNumberOfClusters( new int[]{1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8}, 16 ) );
    }
}
