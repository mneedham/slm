import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static java.util.Arrays.asList;
import static org.junit.Assert.assertEquals;

public class SLMIT
{
    @Test
    public void shouldFindCommunitiesForTheKarateClubNetwork() throws IOException
    {
        // given
        String input = "karate_club_network.txt";
        String output = "output.txt";
        String[] args = new String[]{input, output, "1", "1.0", "3", "1", "10", "0", "1"};

        // When
        ModularityOptimizer.main( args );

        // then
        List<String> communities = readCommunities( output );

        assertEquals(
                asList( "1", "1", "1", "1", "3", "3", "3", "1", "0", "0", "3", "1", "1", "1", "0", "0", "3", "1", "0",
                        "1", "0", "1", "0", "2", "2", "2", "0", "2", "2", "0", "0", "2", "0", "0" ),
                communities );
    }

    private List<String> readCommunities( String output ) throws IOException
    {
        List<String> communities = new ArrayList<>();
        try ( BufferedReader bufferedReader = new BufferedReader( new FileReader( output ) ) )
        {
            String line;
            while ( (line = bufferedReader.readLine()) != null )
            {
                communities.add( line );
            }
        }
        return communities;
    }
}
