import React from 'react';
import Header from './Header';
import './About.css';
import { useNavigate } from 'react-router-dom';

function About() {
  const navigate = useNavigate();
  return (
    <>
      <Header />
      <div className="header-spacer"></div>
      <div className="about-container">
        <main className="content">
          
            <section id = "heading">
            <h1>SUBGRAPH COVERAGE PATTERNS</h1>
            </section>
            <section id="desc">
            <h2>Description</h2>
            <p>
            Subgraph Coverage Patterns (SCP) is a graph mining technique used to uncover frequent, repeating substructures across a dataset made up of many individual graphs. These "graphs" are networks of nodes (e.g., drugs, proteins, or molecules) connected by edges (e.g., interactions or relationships). SCP is designed to systematically search through these graphs to find subgraphs — smaller groups of connected nodes and edges — that appear consistently across multiple graphs in the dataset.
            </p>
            <p>
            In many scientific domains — especially in biomedical research and drug discovery — patterns of interactions often hold more value than isolated data points. By identifying these repetitive patterns, SCP allows researchers to:
            </p>
            <ul>
                <li>
                Detect biologically significant motifs that appear in multiple drug or disease interaction networks.
                </li>
                <li>
                Compare complex drug combinations by analyzing the interaction structure rather than just the list of components.
                </li>
                <li>
                Highlight shared mechanisms between treatments for different conditions, revealing opportunities for drug repurposing.
                </li>
            </ul>
            <p>
            Let’s say you have multiple graphs where each graph represents how a combination of drugs interacts within a biological system. These graphs might vary based on disease, patient group, or experimental condition. SCP helps by:
            </p>
            <ul>
                <li>
                Finding shared substructures in successful drug treatments, which may hint at core components of therapeutic effectiveness.
                </li>
                <li>
                Clustering subgraphs that occur more frequently in cases of positive patient response vs. adverse outcomes.
                </li>
                <li>
                Identifying backbone patterns that can guide the design of new drug combinations, reducing time and trial in experimental phases.
                </li>
            </ul>
            <p>
            Imagine each graph as a map of a city showing streets (edges) and locations (nodes). SCP works like a pattern detector that finds neighborhoods (subgraphs) that keep showing up across different cities. If several cities (graphs) have the same type of neighborhood arrangement — say, a hospital next to a park and a school — then that recurring pattern might be worth studying, especially if it’s associated with better health or outcomes.
            </p>
          </section>

          <section id="param">
            <h2>Key parameters</h2>
            <ul>
              <li>
              Minimum Support Threshold: Ensures that the pattern appears in a sufficient proportion of graphs in the dataset, filtering out infrequent patterns.
              </li>
              <li>
              Maximum Allowed Overlap: Limits how much individual subgraphs within a pattern can share common nodes or edges, preserving pattern distinctiveness.
              </li>
              <li>
              Minimum Coverage Support: Requires the pattern to collectively cover a minimum number of graphs, ensuring the pattern's broad relevance across the dataset.
              </li>
            </ul>
            Guidelines for users - Start off with high minimum supports and relative frequency, then gradually reduce them. For maximum overlap start of low then increase gradually.
          </section>

          <section id="cases">
            <h2>Use cases</h2>
            <ul>
              <li>Drug Repurposing: Identifying recurring drug interaction patterns that work across diseases.</li>
              <li>Adverse Reaction Clustering: Discovering common substructures linked to adverse drug interactions.</li>
              <li>Biomarker Discovery: Spotting graph motifs that frequently occur in disease-specific networks.</li>
              <li>Target Validation: Validating biological pathways by comparing subgraphs across known mechanisms.</li>
            </ul>
            
          </section>


        </main>

        <div className="toc-hover-area">
        <div className="toc-icon">📑</div>
        <div className="hover-toc">
          <h3>Contents</h3>
          <ul>
            <li>
              <a href="#heading">SCP</a>
              <ul>
                <li><a href="#desc">Description</a></li>
                <li><a href="#param">Key Parameters</a></li>
                <li><a href="#cases">Use Cases</a></li>
                <li><a href="#usage">How to Use</a></li>
              </ul>
            </li>
          </ul>
        </div>
      </div>
      </div>

      <button className="try-tool-btn" onClick={() => navigate('/tool')}>Try our Tool</button>
    </>
  );
}

export default About;
