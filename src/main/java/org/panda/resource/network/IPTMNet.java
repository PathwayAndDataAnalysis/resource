package org.panda.resource.network;

import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.FileUtil;
import org.panda.utility.TermCounter;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.SiteSpecificGraph;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Serves the iPTMNet database.
 *
 * To update, go to http://annotation.dbi.udel.edu/iptmnet_data/ and use the file "ptm.txt". Rename it to iPTMNet.txt.
 *
 * @author Ozgun Babur
 */
public class IPTMNet extends FileServer
{
	static IPTMNet instance;

	Map<SignedType, SiteSpecificGraph> graphs;

	public static synchronized IPTMNet get()
	{
		if (instance == null) instance = new IPTMNet();
		return instance;
	}

	public SiteSpecificGraph getGraph(SignedType type)
	{
		return graphs.get(type);
	}

	public Map<SignedType, DirectedGraph> getAllGraphs()
	{
		return new HashMap<>(graphs);
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"iPTMNet.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		graphs = new HashMap<>();
		graphs.put(SignedType.PHOSPHORYLATES, new SiteSpecificGraph("iPTMNet Phosphorylation", SignedType.PHOSPHORYLATES.getTag()));
		graphs.put(SignedType.ACETYLATES, new SiteSpecificGraph("iPTMNet Acetylation", SignedType.ACETYLATES.getTag()));
		graphs.put(SignedType.METHYLATES, new SiteSpecificGraph("iPTMNet Methylation", SignedType.METHYLATES.getTag()));

//		Set<String> set = Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).map(l -> l.split("\t")[0]).collect(Collectors.toSet());
//		System.out.println("set = " + set);

		Map<String, SignedType> modMap = new HashMap<>();
		modMap.put("PHOSPHORYLATION", SignedType.PHOSPHORYLATES);
		modMap.put("ACETYLATION", SignedType.ACETYLATES);
		modMap.put("METHYLATION", SignedType.METHYLATES);

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).map(l -> l.split("\t"))
			.filter(t -> t.length > 7)
			.filter(t -> modMap.keySet().contains(t[0]))
			.filter(t -> t[4].equals("Homo sapiens (Human)"))
			.filter(t -> !t[7].isEmpty())
			.forEach(t ->
		{
			String source = t[7];
			String target = t[3];
			String site = t[5];

			Set<String> pmids = t.length > 9 ? Arrays.stream(t[9].split(",")).map(s -> "PMID:" + s).collect(Collectors.toSet()) : Collections.emptySet();

			SiteSpecificGraph graph = graphs.get(modMap.get(t[0]));

			if (!graph.hasRelation(source, target))
			{
				graph.putRelation(source, target, pmids, site);
			}
			else
			{
				graph.addSite(source, target, site);
				graph.addMediators(source, target, pmids);
			}
		});

		return true;
	}

	public static void main(String[] args)
	{
		BufferedWriter writer = FileUtil.newBufferedWriter("/home/ozgunbabur/Data/IPTMNet.sif");
		for (SignedType type : SignedType.values())
		{
			SiteSpecificGraph graph = get().getGraph(type);
			if (graph != null)
			{
				graph.write(writer);
				graph.printStats();
			}
		}
		FileUtil.closeWriter(writer);
	}
}
