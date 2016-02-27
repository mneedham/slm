LOAD CSV FROM "file:///Users/markneedham/projects/slm/karate_club_network.txt" as row
FIELDTERMINATOR "\t"
MERGE (person1:Person {id: row[0]})
MERGE (person2:Person {id: row[1]})
MERGE (person1)-[:KNOWS]->(person2)