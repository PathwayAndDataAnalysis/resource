package org.panda.resource.network;

import org.biopax.paxtools.pattern.miner.SIFType;
import org.panda.resource.ResourceDirectory;
import org.panda.resource.SignedInteractionText;
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

		Map<String, Map<String, SignedInteractionText>> mapmap = new HashMap<>();

		Files.createDirectories(Paths.get(getPrivateDirectory()));

		// interaction key to sites
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

					SignedInteractionText sit = new SignedInteractionText(line);
					falseMap.put(sit.key(), sit.getSites() == null ? Collections.emptySet() : sit.getSites());
				}
			}
			else
			{
				while (sc.hasNextLine())
				{
					String line = sc.nextLine();
					if (line.isEmpty() || line.startsWith("#")) continue;

					SignedInteractionText sit = new SignedInteractionText(line);

					// Do not consider the relation if it is in the false set
					Set<String> falseSites = falseMap.get(sit.key());
					if (falseSites != null)
					{
						if (falseSites.isEmpty()) continue;
						if (sit.getSites() == null || sit.getSites().isEmpty()) continue;

						sit.getSites().removeAll(falseSites);
						if (sit.getSites().isEmpty()) continue;
					}

					if (!mapmap.containsKey(sit.getType().getTag()))
					{
						mapmap.put(sit.getType().getTag(), new HashMap<>());
					}

					if (mapmap.get(sit.getType().getTag()).containsKey(sit.key()))
					{
						mapmap.get(sit.getType().getTag()).get(sit.key()).merge(sit);
					}
					else
					{
						mapmap.get(sit.getType().getTag()).put(sit.key(), sit);
					}
				}
			}
		}

		for (String type : mapmap.keySet())
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(getPrivateDirectory() + type + ".txt"));

			for (String key : mapmap.get(type).keySet())
			{
				SignedInteractionText sit = mapmap.get(type).get(key);
				writer.write(sit.toStringWOType() + "\n");
			}

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

//		DirectedGraph graph = SignedPCNoTransfac.get().getGraph(SignedType.UPREGULATES_EXPRESSION);
//		Set<String> set = graph.getDownstream("RB1");
//		System.out.println("set.contains(\"MYC\") = " + set.contains("MYC"));
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
