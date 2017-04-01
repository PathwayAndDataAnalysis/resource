package org.panda.resource.network;

import org.biopax.paxtools.pattern.miner.SIFType;
import org.panda.resource.ResourceDirectory;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.CollectionUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.PhosphoGraph;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Serves the signed and directed version of Pathway Commons that contains edges: phosphorylates, dephosphorylates,
 * upregulates-expression, and downregulates-expression.
 *
 * @author Ozgun Babur
 */
public class SignedPC extends PathwayCommons
{
	private static SignedPC instance;

	public static SignedPC get()
	{
		if (instance == null) instance = new SignedPC();
		return instance;
	}

	public Map<SignedType, DirectedGraph> getAllGraphs()
	{
		Map<SignedType, DirectedGraph> map = new HashMap<>();
		for (SignedType type : SignedType.values())
		{
			map.put(type, getGraph(type));
		}
		return map;
	}

	protected String getPrivateDirectory()
	{
		return ResourceDirectory.get() + File.separator + "SignedPC/";
	}

	public DirectedGraph getGraph(SIFType... types)
	{try{
		String edgeType = CollectionUtil.merge(
			Arrays.stream(types).map(SIFType::getTag).collect(Collectors.toList()), ",");

		boolean phos = Arrays.stream(types).anyMatch(type ->
			type == SignedType.PHOSPHORYLATES || type == SignedType.DEPHOSPHORYLATES);

		if (fileExists(types))
		{
			DirectedGraph graph = phos ? new PhosphoGraph("Signed PC", edgeType) :
				new DirectedGraph("Signed PC", edgeType);

			for (SIFType type : types)
			{
				Files.lines(Paths.get(getPrivateDirectory() + type.getTag() + ".txt"))
					.map(line -> line.split("\t")).forEach(token -> {
					if (token.length > 2)
					{
						if (phos && token.length > 3)
						{
							((PhosphoGraph) graph).putRelation(token[0], token[1], token[2], token[3]);
						} else
						{
							graph.putRelation(token[0], token[1], token[2]);
						}
					} else
					{
						graph.putRelation(token[0], token[1]);
					}
				});
			}

			return graph;
		}
		return null;
	} catch (IOException e){throw new RuntimeException(e);}}

	public Graph getSingleGraph(SIFType type)
	{
		return getGraph(type);
	}

	@Override
	public boolean processTheDownloadedFiles()
	{try{
		Map<String, Writer> writers = new HashMap<>();
		Files.createDirectories(Paths.get(getPrivateDirectory()));

		Map<String, Set<String>> falseMap = new HashMap<>();

		for (String localFilename : getLocalFilenames())
		{
			Scanner sc = new Scanner(new File(locateInBase(localFilename)));

			if (localFilename.contains("false"))
			{
				while (sc.hasNextLine())
				{
					String line = sc.nextLine();
					if (line.startsWith("#") || line.isEmpty()) continue;
					String[] token = line.split("\t");
					if (token.length == 3) falseMap.put(line, Collections.emptySet());
					else
					{
						falseMap.put(token[0] + "\t" + token[1] + "\t" + token[2],
							new HashSet<>(Arrays.asList(token[3].split(";"))));
					}
				}
			}
			else
			{
				while (sc.hasNextLine())
				{
					String line = sc.nextLine();
					String[] token = line.split("\t");

					if (token.length > 2)
					{
						// Do not consider the relation if it is in the false set
						String key = token[0] + "\t" + token[1] + "\t" + token[2];
						Set<String> falseSites = falseMap.get(key);
						Set<String> sites = null;
						if (falseSites != null)
						{
							if (falseSites.isEmpty()) continue;
							if (token.length <= 4) continue;
							sites = new HashSet<>(Arrays.asList(token[4].split(";")));
							sites.removeAll(falseSites);
							if (sites.isEmpty()) continue;
						}

						if (!writers.containsKey(token[1])) writers.put(token[1],
							new BufferedWriter(new FileWriter(getPrivateDirectory() + token[1] + ".txt")));

						writers.get(token[1]).write(token[0] + "\t" + token[2]);

						if (token.length > 3)
						{
							writers.get(token[1]).write("\t" + token[3]);
						}
						if (token.length > 4)
						{
							if (falseMap.containsKey(key))
								writers.get(token[1]).write("\t" + CollectionUtil.merge(sites, ";"));
							else writers.get(token[1]).write("\t" + token[4]);
						}

						writers.get(token[1]).write("\n");
					}
				}
			}
		}

		for (Writer writer : writers.values())
		{
			writer.close();
		}

		return true;
	}catch (IOException e){throw new RuntimeException(e);}}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"curated-signed-false.sif", "SignedPC.sif", "curated-signed.sif"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{GITHUB_REPO_BASE + "curated-signed-false.sif", GITHUB_REPO_BASE + "SignedPC.sif.gz",
			GITHUB_REPO_BASE + "curated-signed.sif"};
	}

	public static void main(String[] args)
	{
		printNetworkSizes();
		DirectedGraph graph = SignedPC.get().getGraph(SignedType.PHOSPHORYLATES);
		Set<String> set = graph.getDownstream("TP53");
	}

	private static void printNetworkSizes()
	{
		SignedPC spc = new SignedPC();
		for (SignedType type : SignedType.values())
		{
			spc.getGraph(type).printStats();
		}
	}
}
