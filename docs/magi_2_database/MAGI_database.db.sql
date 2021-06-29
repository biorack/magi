BEGIN TRANSACTION;
CREATE TABLE IF NOT EXISTS "Retro_rules_substrates" (
	"substrate_ID"	INTEGER,
	"retro_rules_ID"	TEXT,
	"retro_rules_smiles"	TEXT,
	"canonical_smiles"	TEXT,
	"inchi_key"	TEXT,
	PRIMARY KEY("substrate_ID")
);
CREATE TABLE IF NOT EXISTS "Retro_rules_to_rhea_reactions" (
	"retro_rules_ID"	TEXT,
	"rhea_ID"	TEXT
);
CREATE TABLE IF NOT EXISTS "Retro_rules_reactions" (
	"reaction_ID"	INTEGER,
	"retro_rules_ID"	TEXT,
	"retro_rules_smarts"	TEXT,
	"substrate_ID"	INTEGER,
	"diameter"	INTEGER,
	PRIMARY KEY("reaction_ID"),
	FOREIGN KEY("substrate_ID") REFERENCES "Retro_rules_substrates"("substrate_ID")
);
CREATE TABLE IF NOT EXISTS "Precomputed_reactions" (
	"precomputed_ID"	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	"molecule_ID"	INTEGER NOT NULL,
	"reaction_ID"	INTEGER NOT NULL,
	"similarity"	NUMERIC NOT NULL,
	FOREIGN KEY("molecule_ID") REFERENCES "Precomputed_molecules"("molecule_ID"),
	FOREIGN KEY("reaction_ID") REFERENCES "Retro_rules_reactions"("reaction_ID")
);
CREATE TABLE IF NOT EXISTS "Precomputed_molecules" (
	"molecule_ID"	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	"smiles"	TEXT NOT NULL,
	"canonical_smiles"	INTEGER NOT NULL,
	"inchi_key"	TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS "Proteins" (
	"rhea_ID"	INTEGER,
	"uniprot_ID"	TEXT,
	"header"	TEXT,
	"protein_ID"	TEXT,
	PRIMARY KEY("uniprot_ID","rhea_ID")
);
COMMIT;
