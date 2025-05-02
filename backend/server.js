
const express = require('express');
const mongoose = require('mongoose');
const cors = require('cors');
const multer = require('multer');
const path = require('path');
const fs = require('fs');
require('dotenv').config();

const app = express();
const PORT = process.env.PORT || 8000;

// Middleware
app.use(cors());
app.use(express.json());

// Connect to MongoDB
mongoose.connect(process.env.MONGODB_URI, {
  useNewUrlParser: true,
  useUnifiedTopology: true,
  dbName: 'graph_algorithm_tool'
})
.then(() => console.log('MongoDB connected successfully'))
.catch(err => console.error('MongoDB connection error:', err));

// Create Schema and Model for uploaded files
const fileSchema = new mongoose.Schema({
  filename: String,
  content: String,
  uploadDate: {
    type: Date,
    default: Date.now
  }
});

const File = mongoose.model('File', fileSchema);

// Schema for algorithm selections
const algorithmSchema = new mongoose.Schema({
  name: String,
  selectedAt: {
    type: Date,
    default: Date.now
  },
  sessionId: String
});

const Algorithm = mongoose.model('Algorithm', algorithmSchema);

// Schema for parameters
const parameterSchema = new mongoose.Schema({
  min_support: Number,
  max_overlap: Number,
  min_coverage: Number,
  sessionId: String,
  setAt: {
    type: Date,
    default: Date.now
  }
});

const Parameter = mongoose.model('Parameter', parameterSchema);

// Setup file storage
const storage = multer.diskStorage({
  destination: function(req, file, cb) {
    const uploadDir = path.join(__dirname, 'uploads');
    if (!fs.existsSync(uploadDir)) {
      fs.mkdirSync(uploadDir, { recursive: true });
    }
    cb(null, uploadDir);
  },
  filename: function(req, file, cb) {
    cb(null, file.originalname);
  }
});

const upload = multer({ 
  storage,
  fileFilter: function(req, file, cb) {
    // Accept only .txt files
    if (file.mimetype !== 'text/plain' && !file.originalname.endsWith('.txt')) {
      return cb(new Error('Only .txt files are accepted'), false);
    }
    cb(null, true);
  }
});

// Global state for the current session (in a production app, you'd use proper session management)
let currentSession = {
  fileId: null,
  algorithm: null,
  parameters: null
};

// Serve static files for datasets
app.use('/datasets', express.static(path.join(__dirname, 'datasets')));

// Routes
app.post('/uploadfile', upload.single('file_upload'), async (req, res) => {
  try {
    if (!req.file) {
      return res.status(400).json({ message: 'No file uploaded' });
    }

    const fileContent = fs.readFileSync(req.file.path, 'utf8');
    
    // Save file to database
    const newFile = new File({
      filename: req.file.originalname,
      content: fileContent
    });
    
    const savedFile = await newFile.save();
    currentSession.fileId = savedFile._id;
    
    res.status(200).json({ 
      message: 'File uploaded successfully',
      fileId: savedFile._id
    });
  } catch (error) {
    console.error('Error uploading file:', error);
    res.status(500).json({ message: 'Error uploading file', error: error.message });
  }
});

app.post('/set-algorithm', async (req, res) => {
  try {
    const { algorithm } = req.body;
    if (!algorithm) {
      return res.status(400).json({ message: 'Algorithm is required' });
    }

    // Save algorithm to database
    const newAlgorithm = new Algorithm({
      name: algorithm,
      sessionId: 'session-' + Date.now() // Simple session ID for demo
    });
    
    const savedAlgorithm = await newAlgorithm.save();
    currentSession.algorithm = savedAlgorithm._id;
    
    res.status(200).json({ 
      message: 'Algorithm set successfully',
      algorithmId: savedAlgorithm._id 
    });
  } catch (error) {
    console.error('Error setting algorithm:', error);
    res.status(500).json({ message: 'Error setting algorithm', error: error.message });
  }
});

app.post('/set-parameters', async (req, res) => {
  try {
    const { min_support, max_overlap, min_coverage } = req.body;
    
    // Validate parameters
    if (min_support === undefined || max_overlap === undefined || min_coverage === undefined) {
      return res.status(400).json({ message: 'All parameters are required' });
    }

    // Save parameters to database
    const newParameters = new Parameter({
      min_support,
      max_overlap,
      min_coverage,
      sessionId: 'session-' + Date.now() // Simple session ID for demo
    });
    
    const savedParameters = await newParameters.save();
    currentSession.parameters = savedParameters._id;
    
    res.status(200).json({ 
      message: 'Parameters set successfully',
      parametersId: savedParameters._id 
    });
  } catch (error) {
    console.error('Error setting parameters:', error);
    res.status(500).json({ message: 'Error setting parameters', error: error.message });
  }
});

app.post('/run', async (req, res) => {
  try {
    // Check if all required data is available
    if (!currentSession.fileId || !currentSession.algorithm || !currentSession.parameters) {
      return res.status(400).json({ 
        message: 'Missing required data to run algorithm',
        missing: {
          file: !currentSession.fileId,
          algorithm: !currentSession.algorithm,
          parameters: !currentSession.parameters
        }
      });
    }

    // Fetch all the data needed to run the algorithm
    const file = await File.findById(currentSession.fileId);
    const algorithm = await Algorithm.findById(currentSession.algorithm);
    const parameters = await Parameter.findById(currentSession.parameters);

    if (!file || !algorithm || !parameters) {
      return res.status(404).json({ message: 'One or more required resources not found' });
    }

    // Here you would implement the actual algorithm based on the parameters
    // For now, we'll return a mock result
    const result = {
      status: 'success',
      algorithm: algorithm.name,
      parameters: {
        min_support: parameters.min_support,
        max_overlap: parameters.max_overlap,
        min_coverage: parameters.min_coverage
      },
      executionTime: Math.random() * 5 + 1, // Random execution time between 1-6 seconds
      coveredGraphs: Math.floor(Math.random() * 10) + 1, // Random number of covered graphs
      discoveredPatterns: []
    };

    // Generate some mock discovered patterns
    for (let i = 0; i < 5; i++) {
      result.discoveredPatterns.push({
        id: `pattern-${i + 1}`,
        support: Math.random().toFixed(2),
        nodes: Math.floor(Math.random() * 5) + 3,
        edges: Math.floor(Math.random() * 8) + 2
      });
    }

    res.status(200).json(result);
  } catch (error) {
    console.error('Error running algorithm:', error);
    res.status(500).json({ message: 'Error running algorithm', error: error.message });
  }
});

// Sample datasets endpoint to list available datasets
app.get('/available-datasets', (req, res) => {
  const datasetsDir = path.join(__dirname, 'datasets');
  
  try {
    if (!fs.existsSync(datasetsDir)) {
      fs.mkdirSync(datasetsDir, { recursive: true });
      // Create sample dataset files for testing if they don't exist
      const sampleDatasets = [
        'graph_transactional_dataset1.txt',
        'graph_transactional_dataset2.txt',
        'graph_transactional_dataset3.txt',
        'graph_transactional_dataset4.txt'
      ];
      
      sampleDatasets.forEach((filename, index) => {
        fs.writeFileSync(
          path.join(datasetsDir, filename),
          `Sample content for dataset ${index + 1}\nThis is a test dataset file.`
        );
      });
    }
    
    const files = fs.readdirSync(datasetsDir)
      .filter(file => file.endsWith('.txt'));
    
    res.status(200).json({ datasets: files });
  } catch (error) {
    console.error('Error listing datasets:', error);
    res.status(500).json({ message: 'Error listing datasets', error: error.message });
  }
});

app.listen(PORT, () => {
  console.log(`Server running on port ${PORT}`);
});

