import React from 'react';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import Home from './Home';
import Citation from './Citation';
import About from './About';
import Tool from './Tool';
import ScrollToTop from './ScrollToTop'; // import the component

function App() {
  return (
    <Router>
      <ScrollToTop /> {/* add this line inside Router */}
      <Routes>
        <Route path="/" element={<Home />} />
        <Route path="/citations" element={<Citation />} />
        <Route path="/about" element={<About />} />
        <Route path="/tool" element={<Tool />} />
      </Routes>
    </Router>
  );
}

export default App;
