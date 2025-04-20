import React, { useState, useRef } from "react";
import "./Tool.css";
import Header from './Header';
import ReactModal from 'react-modal';


const Tool = () => {
  const [step, setStep] = useState(1);
  const [selectedFile, setSelectedFile] = useState(null);
  const [dataset, setDataset] = useState(""); // empty by default
  const [minSupport, setMinSupport] = useState(0.5);
  const [maxOverlap, setMaxOverlap] = useState(0.5);
  const [minCoverage, setMinCoverage] = useState(0.5);
  const [fileText, setFileText] = useState("");
  const [algorithm, setAlgorithm] = useState("");
  const [fetchedText, setFetchedText] = useState("");  // stores content in background
  const [error, setError] = useState("");


  const dropRef = useRef(null);

  const handleNext = () => {
    setStep(prev => prev + 1);
  };

  const handleFileChange = (e) => {
    if (step > 1) return;
    const file = e.target.files[0];
    setSelectedFile(file);
  };

  const handleDrop = (e) => {
    e.preventDefault();
    const file = e.dataTransfer.files[0];
  
    if (!file) return;
  
    if (file.type !== "text/plain" && !file.name.endsWith(".txt")) {
      setError("Only .txt files are accepted.");
      return;
    }
  
    setError("");
    setSelectedFile(file);
    setFileText("");
    const reader = new FileReader();
    reader.onload = (event) => {
      setFileText(event.target.result);
      setFetchedText("");
    };
    reader.readAsText(file);
  };

  const handleRun = async () => {
    let fileContent = "";
  
    if (selectedFile) {
      if (fileText) {
        fileContent = fileText;  // already loaded (drag/drop)
      } else if (fetchedText) {
        fileContent = fetchedText;  // from predefined dataset
      } else {
        const reader = new FileReader();
        reader.onload = async (event) => {
          fileContent = event.target.result;
          await sendToBackend(fileContent);
        };
        reader.readAsText(selectedFile);
        return;  // exit to wait for reader callback
      }
    }
  
    await sendToBackend(fileContent);
  };
  
  const sendToBackend = async (content) => {
    try {
      const payload = {
        algorithm,
        min_support: minSupport,
        max_overlap: maxOverlap,
        min_coverage: minCoverage,
        file_content: content,
        filename: selectedFile?.name || dataset,
      };
  
      const response = await fetch('/api/run-algorithm', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(payload),
      });
  
      if (!response.ok) {
        throw new Error(`Server error: ${response.statusText}`);
      }
  
      const result = await response.json();
      alert("Algorithm completed successfully.");
      console.log(result); // You can display results in the UI later.
    } catch (err) {
      console.error("Error sending data to backend:", err);
      alert("Error: " + err.message);
    }
  };
  


  const [showModal, setShowModal] = useState(false);

  const closeModal = () => setShowModal(false);

  const handleViewFile = () => {
    if (!selectedFile) return;

    if (selectedFile.name.startsWith("graph_transactional_dataset")) {
      // Preset file - preview from fetchedText
      setFileText(fetchedText);
      setShowModal(true);
    } else {
      // Uploaded file - use FileReader
      const reader = new FileReader();
      reader.onload = (e) => {
        setFileText(e.target.result);
        setShowModal(true);
      };
      reader.onerror = (e) => alert("Error reading file: " + e.target.error.message);
      reader.readAsText(selectedFile);
    }
  };

  


  return (
    <>
      <Header />
      <div className="tool-container">

        {/* Step 1 */}
        <h2>1 Upload datasets</h2>
        <div
          className={`upload-box ${step > 1 ? "disabled" : ""}`}
          ref={dropRef}
          onDrop={handleDrop}
          onDragOver={(e) => e.preventDefault()}
        >
          <label htmlFor="file-upload" className="file-upload-label">
            <div className="upload-icon">☁️</div>
            <div>
              {selectedFile ? (
                <div className="file-display">
                  <span>{selectedFile.name}</span>
                  <button
                    className="view-button"
                    onClick={handleViewFile}
                    disabled={step > 1}
                  >
                    View Dataset
                  </button>
                </div>
              ) : (
                "Drag and drop a txt file here, or click to select a txt file"
              )}
            </div>
          </label>
          <input
            id="file-upload"
            type="file"
            accept=".txt"
            onChange={handleFileChange}
            className="file-input"
            disabled={step > 1}
          />
          <div className="dropdown-group">
            <span>Or use available ones</span>
            <select
              value={dataset}
              onChange={async (e) => {
                const selected = e.target.value;

                if (selected === "") {
                  setDataset("");
                  setSelectedFile(null);
                  return;
                }

                setDataset(selected);
                setSelectedFile({ name: selected });

                try {
                  const response = await fetch(`/datasets/${selected}`);
                  const text = await response.text();
                  setFetchedText(text);
                  setFileText("");
                } catch (err) {
                  console.error("Error loading dataset:", err);
                  setFetchedText("");
                  setFileText("");
                  setSelectedFile(null);
                }
              }}
              className="dropdown"
              disabled={step > 1}
            >
              <option value="">-- Select a dataset --</option>
              <option value="graph_transactional_dataset1.txt">graph_transactional_dataset1.txt</option>
              <option value="graph_transactional_dataset2.txt">graph_transactional_dataset2.txt</option>
              <option value="graph_transactional_dataset3.txt">graph_transactional_dataset3.txt</option>
              <option value="graph_transactional_dataset4.txt">graph_transactional_dataset4.txt</option>
            </select>
          </div>
        </div>
        {error && <p className="error-message">{error}</p>}

        {/* Only show the Next button if a file is selected */}
        {step === 1 && selectedFile && (
          <button className="next-button" onClick={handleNext}>
            Next
          </button>
        )}


        {/* Step 2 */}
        {step >= 2 && (
          <>
            <h2>2 Select Algorithm</h2>
            <select
              value={algorithm}
              onChange={(e) => setAlgorithm(e.target.value)}
              className="dropdown"
              disabled={step > 2}
            >
              <option value="" disabled>-- Select an algorithm --</option>
              <option value="Subgraph Coverage Patterns">Subgraph Coverage Patterns</option>
            </select>
            {step === 2 && algorithm && <button className="next-button" onClick={handleNext}>Next</button>}
          </>
        )}

        {/* Step 3 */}
        {step >= 3 && (
          <>
            <h2>3 Set parameters</h2>
            <div className="slider-container">
              {/* Min Support */}
              <label>a. Minimum Relative Frequency</label>
              <input type="range" min="0" max="1" step="0.01" value={minSupport}
                onChange={(e) => setMinSupport(Number(e.target.value))} disabled={step > 3} />
              <input type="number" min="0" max="1" step="0.01" value={minSupport}
                onChange={(e) => setMinSupport(Number(e.target.value))} disabled={step > 3} />
              <p className="slider-info">Minimum percentage of graphs a subgraph must appear in to be considered frequent.</p>

              {/* Max Overlap */}
              <label>b. Maximum Allowed Overlap</label>
              <input type="range" min="0" max="1" step="0.01" value={maxOverlap}
                onChange={(e) => setMaxOverlap(Number(e.target.value))} disabled={step > 3} />
              <input type="number" min="0" max="1" step="0.01" value={maxOverlap}
                onChange={(e) => setMaxOverlap(Number(e.target.value))} disabled={step > 3} />
              <p className="slider-info">Maximum percentage of overlap allowed between discovered subgraphs.</p>

              {/* Min Coverage */}
              <label>c. Minimum Coverage Support</label>
              <input type="range" min="0" max="1" step="0.01" value={minCoverage}
                onChange={(e) => setMinCoverage(Number(e.target.value))} disabled={step > 3} />
              <input type="number" min="0" max="1" step="0.01" value={minCoverage}
                onChange={(e) => setMinCoverage(Number(e.target.value))} disabled={step > 3} />
              <p className="slider-info">Minimum number of graphs that the selected subgraphs should collectively cover.</p>
              <p className="guidelines">Guidelines: Start off with high minimum supports and relative frequency, then gradually reduce them. For maximum overlap start of low then increase gradually.</p>
            </div>
            {step === 3 && <button className="next-button" onClick={handleNext}>Next</button>}
          </>
        )}

        {/* Step 4 */}
        {step >= 4 && (
          <>
            <h2>4 View results</h2>
            <button className="run-button" onClick={handleRun}>Run Algorithm</button>
          </>
        )}
      </div>
      <ReactModal
        isOpen={showModal}
        onRequestClose={closeModal}
        contentLabel="Dataset Preview"
        className="Modal"
        overlayClassName="Overlay"
      >
        <div className="modal-header">
          <h3>Dataset Preview: {selectedFile?.name}</h3>
          <button onClick={closeModal} className="close-button">X</button>
        </div>
        <pre className="modal-content">{fileText}</pre>
      </ReactModal>
    </>
  );
};

export default Tool;
