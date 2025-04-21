import React, { useState } from 'react';
import { Link, useLocation } from 'react-router-dom';
import './Header.css';

function Header() {
  const location = useLocation();
  const [menuOpen, setMenuOpen] = useState(false);

  const toggleMenu = () => setMenuOpen(!menuOpen);

  return (
    <header>
      <Link to="/" className="logo">
        <h1>SCPMiner</h1>
      </Link>
      <button className="menu-toggle" onClick={toggleMenu}>
        ☰
      </button>
      <nav className={menuOpen ? 'open' : ''}>
        <Link to="/" className={location.pathname === '/' ? 'active' : ''}>Home</Link>
        <Link to="/about" className={location.pathname === '/about' ? 'active' : ''}>About</Link>
        <Link to="/tool" className={location.pathname === '/tool' ? 'active' : ''}>Tool</Link>
        <Link to="/citations" className={location.pathname === '/citations' ? 'active' : ''}>Citations</Link>
      </nav>
    </header>
  );
}

export default Header;