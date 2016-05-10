package org.panda.resource.network;

import org.biopax.paxtools.pattern.miner.SIFType;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.CollectionUtil;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.PhosphoGraph;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class SignedPC extends PathwayCommons
{
	public Map<SignedType, Graph> getAllGraphs()
	{
		Map<SignedType, Graph> map = new HashMap<>();
		for (SignedType type : SignedType.values())
		{
			map.put(type, getGraph(type));
		}
		return map;
	}

	protected String getPrivateDirectory()
	{
		return getResourceFilesLocation() + File.separator + "SignedPC/";
	}

	public Graph getGraph(SIFType... types)
	{try{
		String edgeType = CollectionUtil.merge(
			Arrays.stream(types).map(SIFType::getTag).collect(Collectors.toList()), ",");

		boolean phos = Arrays.stream(types).anyMatch(type ->
			type == SignedType.PHOSPHORYLATES || type == SignedType.DEPHOSPHORYLATES);

		if (fileExists(types))
		{
			Graph graph = phos ? new PhosphoGraph("Signed PC", edgeType) :
				new Graph("Signed PC", edgeType);

			for (SIFType type : types)
			{
				Files.lines(Paths.get(getPrivateDirectory() + type.getTag() + ".txt"))
					.map(line -> line.split("\t")).forEach(token -> {
					if (token.length > 2)
					{
						if (phos && token.length > 3)
						{
							((PhosphoGraph) graph).putRelation(token[0], token[1], token[2], type.isDirected(), token[3]);
						} else
						{
							graph.putRelation(token[0], token[1], token[2], type.isDirected());
						}
					} else
					{
						graph.putRelation(token[0], token[1], type.isDirected());
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
		Scanner sc = new Scanner(new File(locateInBase(getLocalFilenames()[0])));

		Map<String, Writer> writers = new HashMap<String, Writer>();

		Files.createDirectories(Paths.get(getPrivateDirectory()));

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] token = line.split("\t");

			if (token.length > 2)
			{
				if (!writers.containsKey(token[1])) writers.put(token[1],
					new BufferedWriter(new FileWriter(getPrivateDirectory() + token[1] + ".txt")));

				writers.get(token[1]).write(token[0] + "\t" + token[2]);

				if (token.length > 3)
				{
					writers.get(token[1]).write("\t" + token[3]);
				}
				if (token.length > 4)
				{
					writers.get(token[1]).write("\t" + token[4]);
				}

				writers.get(token[1]).write("\n");
			}
		}

		for (Writer writer : writers.values())
		{
			writer.close();
		}

		return Files.deleteIfExists(Paths.get(locateInBase(getLocalFilenames()[0])));
	}catch (IOException e){throw new RuntimeException(e);}}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"SignedPC.sif"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{GITHUB_REPO_BASE + "SignedPC.sif.gz"};
	}

	public static void main(String[] args)
	{
		printNetworkSizes();
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
