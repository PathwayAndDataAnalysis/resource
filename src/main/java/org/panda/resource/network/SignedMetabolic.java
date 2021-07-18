package org.panda.resource.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.DirectedGraph;

public class SignedMetabolic
{
	public static DirectedGraph getGraph(SignedType type)
	{
		switch (type)
		{
			case PRODUCES:
				DirectedGraph graph1 = (DirectedGraph) PathwayCommons.get().getGraph(SIFEnum.CONTROLS_PRODUCTION_OF);
				graph1.setEdgeType(SignedType.PRODUCES.getTag());
				return graph1;
			case USED_TO_PRODUCE:
				DirectedGraph graph2 = (DirectedGraph) PathwayCommons.get().getGraph(SIFEnum.USED_TO_PRODUCE);
				graph2.setEdgeType(SignedType.USED_TO_PRODUCE.getTag());
				return graph2;
			case CONSUMES:
			{
				DirectedGraph graph3 = (DirectedGraph) PathwayCommons.get().getGraph(SIFEnum.CONSUMPTION_CONTROLLED_BY);
				graph3 = graph3.getReverseGraph();
				graph3.setEdgeType(SignedType.CONSUMES.getTag());
				return graph3;
			}
			default: throw new IllegalArgumentException("The edge type " + type + " is not supported.");
		}
	}

	public static void main(String[] args)
	{
		DirectedGraph graph = getGraph(SignedType.USED_TO_PRODUCE);
		String source = graph.getOneSideSymbols(true).iterator().next();
		String target = graph.getDownstream(source).iterator().next();
		System.out.println("source = " + source);
		System.out.println("target = " + target);
	}
}
