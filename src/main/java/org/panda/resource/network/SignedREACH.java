package org.panda.resource.network;

import org.panda.resource.ResourceDirectory;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Serves the signed and directed version of REACH that contains edges: phosphorylates and dephosphorylates.
 *
 * @author Ozgun Babur
 */
public class SignedREACH extends SignedPC
{
	private static SignedREACH instance;

	public static synchronized SignedREACH get()
	{
		if (instance == null) instance = new SignedREACH();
		return instance;
	}

	public Map<SignedType, DirectedGraph> getAllGraphs()
	{
		Map<SignedType, DirectedGraph> map = new HashMap<>();
		for (SignedType type : SignedType.values())
		{
			if (type.isSiteSpecific())
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
		DirectedGraph graph = SignedREACH.get().getGraph(SignedType.PHOSPHORYLATES);
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
