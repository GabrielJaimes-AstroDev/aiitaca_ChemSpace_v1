import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
import re

# Page configuration
st.set_page_config(
    page_title="Chemical Space Visualization",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #2ca02c;
        margin-bottom: 1rem;
    }
    .info-box {
        background-color: #f0f8ff;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
        margin-bottom: 1rem;
    }
    .molecule-card {
        background-color: #f8f9fa;
        padding: 0.5rem;
        border-radius: 0.3rem;
        margin: 0.2rem 0;
        border-left: 3px solid #9467bd;
    }
</style>
""", unsafe_allow_html=True)

def load_plot_data(file_path):
    """Load plot data from JSON file"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        st.error(f"Error loading plot data: {str(e)}")
        return None

def format_chemical_formula(formula):
    """Format chemical formula with subscripts"""
    if not formula or pd.isna(formula):
        return ""
    
    # Handle subscripts
    formula = re.sub(r'(\d+)', r'<sub>\1</sub>', str(formula))
    
    # Handle ionic charges
    formula = re.sub(r'([+-]\d*)', r'<sup>\1</sup>', formula)
    
    return formula

def create_interactive_plot(plot_data):
    """Create interactive Plotly visualization"""
    
    # Create DataFrame from plot data
    df_points = pd.DataFrame({
        'x': [point[0] for point in plot_data['points']],
        'y': [point[1] for point in plot_data['points']],
        'label': plot_data['labels'],
        'database': plot_data['databases'],
        'formula': plot_data['formulas'],
        'color': plot_data['colors']
    })
    
    # Create the main scatter plot
    fig = go.Figure()
    
    # Add points for each database
    databases = df_points['database'].unique()
    for db in databases:
        db_mask = df_points['database'] == db
        db_data = df_points[db_mask]
        
        # Create hover text
        hover_text = []
        for _, row in db_data.iterrows():
            text = f"Database: {row['database']}<br>"
            if row['label']:
                text += f"Name: {row['label']}<br>"
            if row['formula']:
                formatted_formula = format_chemical_formula(row['formula'])
                text += f"Formula: {formatted_formula}"
            hover_text.append(text)
        
        fig.add_trace(go.Scatter(
            x=db_data['x'],
            y=db_data['y'],
            mode='markers',
            marker=dict(
                color=db_data['color'].iloc[0],
                size=8 if db == 'PubChem' else 12,
                opacity=0.3 if db == 'PubChem' else 0.7
            ),
            name=db,
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=True
        ))
    
    # Add KNN connections if they exist
    if plot_data.get('knn_connections'):
        connections = plot_data['knn_connections']
        for connection in connections:
            guapos_idx = connection['guapos_index']
            neighbor_idx = connection['neighbor_index']
            
            if (guapos_idx < len(plot_data['points']) and 
                neighbor_idx < len(plot_data['points'])):
                
                x0, y0 = plot_data['points'][guapos_idx]
                x1, y1 = plot_data['points'][neighbor_idx]
                
                fig.add_trace(go.Scatter(
                    x=[x0, x1],
                    y=[y0, y1],
                    mode='lines',
                    line=dict(color='gray', width=1, dash='dash'),
                    opacity=0.6,
                    showlegend=False,
                    hoverinfo='skip'
                ))
    
    # Update layout
    fig.update_layout(
        title="Chemical Space Visualization - PCA/UMAP Projection",
        xaxis_title="Dimension 1",
        yaxis_title="Dimension 2",
        hovermode='closest',
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        height=700,
        width=1000
    )
    
    fig.update_xaxes(showgrid=False, zeroline=False, showticklabels=False)
    fig.update_yaxes(showgrid=False, zeroline=False, showticklabels=False)
    
    return fig

def create_database_summary(plot_data):
    """Create database summary statistics"""
    databases = plot_data['databases']
    db_counts = pd.Series(databases).value_counts().reset_index()
    db_counts.columns = ['Database', 'Count']
    
    fig = px.bar(
        db_counts,
        x='Database',
        y='Count',
        color='Database',
        title="Molecule Distribution by Database",
        color_discrete_map={
            'GUAPOS': '#1f77b4',
            'TMC_1': '#ff7f0e',
            'All_Discoveries': '#2ca02c',
            'KIDA': '#FFFF00',
            'PubChem': '#9467bd'
        }
    )
    
    fig.update_layout(
        xaxis_title="Database",
        yaxis_title="Number of Molecules",
        showlegend=False
    )
    
    return fig

def create_molecule_info_panel(plot_data, selected_point):
    """Create molecule information panel"""
    if selected_point is None:
        return None
    
    point_idx = selected_point['pointIndex']
    
    if point_idx >= len(plot_data['points']):
        return None
    
    info = {
        'Database': plot_data['databases'][point_idx],
        'Name': plot_data['labels'][point_idx] or 'Unknown',
        'Formula': plot_data['formulas'][point_idx] or 'Unknown',
        'Coordinates': f"({plot_data['points'][point_idx][0]:.3f}, {plot_data['points'][point_idx][1]:.3f})"
    }
    
    return info

def main():
    st.markdown('<h1 class="main-header">🧪 Chemical Space Visualization</h1>', unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.header("📊 Data Configuration")
        
        # File upload
        uploaded_file = st.file_uploader(
            "Upload plot data JSON",
            type=['json'],
            help="Upload the plot_data.json file generated by the analysis script"
        )
        
        if uploaded_file is not None:
            try:
                plot_data = json.load(uploaded_file)
                st.success("✅ Data loaded successfully!")
            except Exception as e:
                st.error(f"❌ Error loading file: {str(e)}")
                plot_data = None
        else:
            # Try to load from default path
            try:
                plot_data = load_plot_data("plot_data.json")
                if plot_data:
                    st.info("📁 Using default plot_data.json file")
            except:
                plot_data = None
                st.warning("⚠️ Please upload a plot data JSON file")
        
        if plot_data:
            st.markdown("---")
            st.header("🎨 Display Options")
            
            # Show database summary
            db_counts = pd.Series(plot_data['databases']).value_counts()
            st.write("**Database Summary:**")
            for db, count in db_counts.items():
                st.write(f"• {db}: {count} molecules")
            
            # Show KNN info if available
            if plot_data.get('knn_connections'):
                st.write(f"• KNN Connections: {len(plot_data['knn_connections'])}")
            
            st.markdown("---")
            st.header("🔍 Interaction Guide")
            st.info("""
            - **Hover** over points to see molecule details
            - **Click** on points to view detailed information
            - **Zoom** with mouse wheel or touchpad
            - **Pan** by dragging the plot
            - **Reset view** with home button in toolbar
            """)
    
    # Main content
    if plot_data is None:
        st.warning("""
        ## Welcome to Chemical Space Visualization!
        
        To get started:
        
        1. Run the main analysis script to generate `plot_data.json`
        2. Upload the JSON file using the sidebar
        3. Or place `plot_data.json` in the same directory as this app
        
        The visualization will show:
        - Molecular distributions across different databases
        - PCA/UMAP projections of chemical space
        - KNN connections between molecules
        - Interactive exploration capabilities
        """)
        
        # Show example image or placeholder
        st.image("https://via.placeholder.com/800x400/1f77b4/ffffff?text=Chemical+Space+Visualization", 
                use_column_width=True)
        
        return
    
    # Create tabs for different views
    tab1, tab2, tab3 = st.tabs(["📈 Interactive Map", "📊 Statistics", "🔍 Molecule Explorer"])
    
    with tab1:
        st.markdown('<h2 class="sub-header">Interactive Chemical Space Map</h2>', unsafe_allow_html=True)
        
        # Create interactive plot
        fig = create_interactive_plot(plot_data)
        
        # Display the plot
        st.plotly_chart(fig, use_container_width=True)
        
        # Add some context
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            <div class="info-box">
            <h4>About this Visualization</h4>
            <p>This map shows molecules projected into 2D space using PCA and UMAP dimensionality reduction techniques. 
            Colors represent different databases, and connections show molecular similarities detected by KNN analysis.</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown("""
            <div class="info-box">
            <h4>How to Interpret</h4>
            <p>• <strong>Closer points</strong> = More similar molecules<br>
            • <strong>Colors</strong> = Source database<br>
            • <strong>Gray lines</strong> = KNN similarity connections<br>
            • <strong>Larger points</strong> = GUAPOS molecules</p>
            </div>
            """, unsafe_allow_html=True)
    
    with tab2:
        st.markdown('<h2 class="sub-header">Database Statistics</h2>', unsafe_allow_html=True)
        
        # Create statistics visualizations
        col1, col2 = st.columns(2)
        
        with col1:
            # Database distribution
            db_fig = create_database_summary(plot_data)
            st.plotly_chart(db_fig, use_container_width=True)
        
        with col2:
            # Additional statistics
            total_molecules = len(plot_data['points'])
            unique_databases = len(set(plot_data['databases']))
            
            st.metric("Total Molecules", f"{total_molecules:,}")
            st.metric("Unique Databases", unique_databases)
            
            if plot_data.get('knn_connections'):
                st.metric("KNN Connections", len(plot_data['knn_connections']))
            
            # Show database details
            st.write("**Database Details:**")
            db_counts = pd.Series(plot_data['databases']).value_counts()
            for db, count in db_counts.items():
                percentage = (count / total_molecules) * 100
                st.write(f"• {db}: {count} ({percentage:.1f}%)")
    
    with tab3:
        st.markdown('<h2 class="sub-header">Molecule Explorer</h2>', unsafe_allow_html=True)
        
        # Create a searchable table of molecules
        molecules_df = pd.DataFrame({
            'Database': plot_data['databases'],
            'Name': plot_data['labels'],
            'Formula': plot_data['formulas'],
            'X': [p[0] for p in plot_data['points']],
            'Y': [p[1] for p in plot_data['points']]
        })
        
        # Search and filter options
        col1, col2 = st.columns(2)
        
        with col1:
            selected_db = st.multiselect(
                "Filter by Database",
                options=sorted(molecules_df['Database'].unique()),
                default=sorted(molecules_df['Database'].unique())
            )
        
        with col2:
            search_term = st.text_input("Search by Name or Formula", "")
        
        # Apply filters
        filtered_df = molecules_df[molecules_df['Database'].isin(selected_db)]
        if search_term:
            filtered_df = filtered_df[
                filtered_df['Name'].str.contains(search_term, case=False, na=False) |
                filtered_df['Formula'].str.contains(search_term, case=False, na=False)
            ]
        
        # Display results
        st.write(f"**Found {len(filtered_df)} molecules**")
        
        # Pagination
        page_size = 20
        total_pages = max(1, (len(filtered_df) + page_size - 1) // page_size)
        page = st.number_input("Page", min_value=1, max_value=total_pages, value=1)
        
        start_idx = (page - 1) * page_size
        end_idx = min(start_idx + page_size, len(filtered_df))
        
        # Display molecules
        for idx in range(start_idx, end_idx):
            molecule = filtered_df.iloc[idx]
            with st.expander(f"{molecule['Name'] or 'Unknown'} - {molecule['Database']}"):
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write(f"**Database:** {molecule['Database']}")
                    st.write(f"**Name:** {molecule['Name'] or 'Unknown'}")
                    st.write(f"**Formula:** {format_chemical_formula(molecule['Formula'])}", unsafe_allow_html=True)
                
                with col2:
                    st.write(f"**Coordinates:** ({molecule['X']:.3f}, {molecule['Y']:.3f})")
                    st.write(f"**Index:** {filtered_df.index[idx]}")
        
        # Download option
        csv = filtered_df.to_csv(index=False)
        st.download_button(
            label="📥 Download Filtered Data",
            data=csv,
            file_name="filtered_molecules.csv",
            mime="text/csv"
        )

if __name__ == "__main__":
    main()
