import React from 'react';
import './Citation.css';
import Header from './Header';

const Citation = () => {
  return (
    <>
      <Header />
      <section className="citation-section">
        <div className="container">
          <h2 className="citation-heading">Citations</h2>
          <p className="citation-text">
            Reddy, A. S., Reddy, P. K., Mondal, A., & Priyakumar, U. D. (2021). Mining subgraph coverage patterns from graph transactions. 
            <em> International Journal of Data Science and Analytics, 13(2), 105–121.</em>{' '}
            <a 
              href="https://doi.org/10.1007/s41060-021-00292-y" 
              target="_blank" 
              rel="noopener noreferrer"
            >
              https://doi.org/10.1007/s41060-021-00292-y
            </a>
          </p>
        </div>
      </section>
    </>
  );
};

export default Citation;
