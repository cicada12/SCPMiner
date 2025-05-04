import React from 'react';
import './Citation.css';
import Header from './Header';
const Citation = () => {
  return (
    <>
      <Header />
      <section className="citation-section">
        <div className="container">
          <div className='step-box'>
          <h2 className="citation-heading">Citations</h2>
          <ul className="citation-list">

            <li>
              Reddy, A. S., Reddy, P. K., Mondal, A., & Priyakumar, U. D. (2021). A Model of Graph Transactional Coverage Patterns with Applications to Drug Discovery. 
              <em>2021 IEEE 28th International Conference on High Performance Computing, Data, and Analytics (HiPC).</em>
            </li>
            <li>
              Snowden, M., & Green, D. V. (2008). The impact of diversity-based, high-throughput screening on drug discovery: chance favours the prepared mind. 
              <em>Current Opinion in Drug Discovery & Development, 11(4), 553–558.</em>
            </li>
            <li>
              Aida, M., Pieter, M., Wout, B., Pieter, M., Boris, C., & Goethals, K. L. (2018). Grasping frequent subgraph mining for bioinformatics applications. 
              <em>BioData Mining, 11(1), 1–20.</em>
            </li>
          </ul>
        </div>
        </div>
      </section>
    </>
  );
};

export default Citation;
