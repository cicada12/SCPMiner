import React from 'react';
import Header from './Header';
import './Help.css';

const Help = () => {
  return (
    <>
      <Header />
      <section id="usage">
        <h2>How to Use</h2>
        <p>
          Follow along the tutorial to upload your code files and receive detailed results.
        </p>
        
        <div className="video-wrapper">
          <div className="video-container">
            <iframe
              src="https://www.youtube.com/embed/VeaAMwpiLbI"
              title="Tutorial Video"
              allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
              allowFullScreen
            ></iframe>
          </div>
        </div>
      </section>
    </>
  );
};

export default Help;
