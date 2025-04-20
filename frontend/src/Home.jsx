import React from 'react';
import './Home.css';
import Header from './Header';
import { useNavigate } from 'react-router-dom';


export default function SCPMiner() {
  const navigate = useNavigate();
  return (
    <div className="font-sans">
      {/* Fixed Header */}
      <Header />

      {/* First Section - Full screen content below fixed header */}
      <section className="first_section">
        <div className="box1">
          <h2 className="tagline">A tool for mining SCPs from graph data</h2>
          <p className="scp_description">
            A Subgraph Coverage Pattern is a subgraph pattern whose subgraphs collectively meet the minimum relative
            frequency threshold, satisfy the maximum allowed overlap, and achieve the minimum coverage
            support.
            Subgraph Coverage Patterns can be used in drug discovery to identify frequent molecular substructures
            associated with specific biological activities or therapeutic properties.
          </p>
          <button className="get_started_btn" onClick={() => navigate('/about')}>
            Get started
          </button>
        </div>
      </section>


      {/* Scrollable Second Section */}
      <section className="second_section">
        <div className="second_box">
          <div>
            <h3 className="algo1" onClick={() => navigate('/about')}>Subgraph Coverage Patterns</h3>
            <p className="definition">
              Subgraph coverage pattern identifies frequent subgraphs within a dataset,
              capturing recurring structural patterns in graph-based data. It helps reveal
              relationships, behaviors, or anomalies in complex networks.
            </p>
          </div>
        </div>
      </section>
    </div>
  );
}
