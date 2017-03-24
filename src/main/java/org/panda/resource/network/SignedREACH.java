package org.panda.resource.network;

import org.biopax.paxtools.pattern.miner.SIFType;
import org.panda.resource.ResourceDirectory;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.CollectionUtil;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.PhosphoGraph;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Serves the signed and directed version of REACH that contains edges: phosphorylates and dephosphorylates.
 *
 * @author Ozgun Babur
 */
public class SignedREACH extends SignedPC
{
	private static SignedREACH instance;

	public static SignedREACH get()
	{
		if (instance == null) instance = new SignedREACH();
		return instance;
	}

	public Map<SignedType, Graph> getAllGraphs()
	{
		Map<SignedType, Graph> map = new HashMap<>();
		for (SignedType type : SignedType.values())
		{
			if (type.isPhospho())
			{
				map.put(type, getGraph(type));
			}
		}
		return map;
	}

	protected String getPrivateDirectory()
	{
		return ResourceDirectory.get() + File.separator + "SignedREACH/";
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"SignedREACH.sif"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{GITHUB_REPO_BASE + "SignedREACH.sif.gz"};
	}

	public static void main(String[] args)
	{
		printNetworkSizes();
		Graph graph = SignedREACH.get().getGraph(SignedType.PHOSPHORYLATES);
		Set<String> set = graph.getDownstream("TP53");
	}

	private static void printNetworkSizes()
	{
		SignedREACH spc = new SignedREACH();
		for (SignedType type : SignedType.values())
		{

			spc.getGraph(type).printStats();
		}
	}
}
