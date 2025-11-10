#include "writers/hdf5_writer.hpp"
#include <cstring>
#include <memory>

namespace sph {

HDF5Writer::HDF5Writer(int compression)
    : m_file_id(-1)
    , m_compression(compression)
{
    // Clamp compression level to valid range
    if (m_compression < 0) m_compression = 0;
    if (m_compression > 9) m_compression = 9;
}

HDF5Writer::~HDF5Writer() {
    close();
}

bool HDF5Writer::open(const std::string& filepath, const OutputMetadata& metadata) {
    // Close existing file if open
    if (m_file_id >= 0) {
        close();
    }
    
    // Create new HDF5 file (overwrite if exists)
    m_file_id = H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (m_file_id < 0) {
        return false;
    }
    
    // Write metadata group
    if (!write_metadata_group(metadata)) {
        close();
        return false;
    }
    
    // Write energy group
    if (!write_energy_group(metadata)) {
        close();
        return false;
    }
    
    return true;
}

bool HDF5Writer::write_particles(const std::vector<SPHParticle*>& particles) {
    if (m_file_id < 0) {
        return false;
    }
    
    return write_particles_group(particles);
}

void HDF5Writer::close() {
    if (m_file_id >= 0) {
        H5Fclose(m_file_id);
        m_file_id = -1;
    }
}

bool HDF5Writer::write_metadata_group(const OutputMetadata& metadata) {
    // Create metadata group
    hid_t group_id = H5Gcreate(m_file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0) return false;
    
    // Write attributes as scalars or strings
    hsize_t dims[1] = {1};
    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    
    // Write simple scalar attributes
    hid_t attr_id;
    
    // format_version
    attr_id = H5Acreate(group_id, "format_version", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, &metadata.format_version);
    H5Aclose(attr_id);
    
    // step
    attr_id = H5Acreate(group_id, "step", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, &metadata.step);
    H5Aclose(attr_id);
    
    // particle_count
    attr_id = H5Acreate(group_id, "particle_count", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, &metadata.particle_count);
    H5Aclose(attr_id);
    
    // time_code
    attr_id = H5Acreate(group_id, "time_code", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &metadata.time_code);
    H5Aclose(attr_id);
    
    // time_physical
    attr_id = H5Acreate(group_id, "time_physical", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &metadata.time_physical);
    H5Aclose(attr_id);
    
    // gamma
    attr_id = H5Acreate(group_id, "gamma", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &metadata.gamma);
    H5Aclose(attr_id);
    
    H5Sclose(dataspace_id);
    
    // Write string attributes
    write_string_attribute(group_id, "simulation_name", metadata.simulation_name);
    write_string_attribute(group_id, "timestamp", metadata.timestamp);
    write_string_attribute(group_id, "sph_type", metadata.get_sph_type_name());
    write_string_attribute(group_id, "kernel_type", metadata.get_kernel_type_name());
    
    // Write unit system as JSON string
    nlohmann::json units_json = metadata.units.to_json();
    write_string_attribute(group_id, "unit_system", units_json.dump(2));
    
    // Write full metadata as JSON string for completeness
    nlohmann::json full_meta = metadata.to_json();
    write_string_attribute(group_id, "metadata_json", full_meta.dump(2));
    
    H5Gclose(group_id);
    return true;
}

bool HDF5Writer::write_particles_group(const std::vector<SPHParticle*>& particles) {
    const hsize_t N = particles.size();
    if (N == 0) return false;
    
    // Create particles group
    hid_t group_id = H5Gcreate(m_file_id, "/particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0) return false;
    
    // Setup compression
    hid_t plist_id = H5P_DEFAULT;
    if (m_compression > 0) {
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {N < 1024 ? N : 1024};  // Chunk size
        H5Pset_chunk(plist_id, 1, chunk_dims);
        H5Pset_deflate(plist_id, m_compression);
    }
    
    // Helper lambda to write 1D dataset
    auto write_1d_dataset = [&](const char* name, const std::vector<real>& data) -> bool {
        hsize_t dims[1] = {N};
        hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
        hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_DOUBLE, dataspace_id,
                                     H5P_DEFAULT, plist_id, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        return true;
    };
    
    // Helper lambda to write 1D integer dataset
    auto write_1d_int_dataset = [&](const char* name, const std::vector<int>& data) -> bool {
        hsize_t dims[1] = {N};
        hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
        hid_t plist_int = H5P_DEFAULT;
        if (m_compression > 0) {
            plist_int = H5Pcreate(H5P_DATASET_CREATE);
            hsize_t chunk_dims[1] = {N < 1024 ? N : 1024};
            H5Pset_chunk(plist_int, 1, chunk_dims);
            H5Pset_deflate(plist_int, m_compression);
        }
        hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_INT, dataspace_id,
                                     H5P_DEFAULT, plist_int, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(dataset_id);
        if (plist_int != H5P_DEFAULT) H5Pclose(plist_int);
        H5Sclose(dataspace_id);
        return true;
    };
    
    // Extract data into vectors
    std::vector<int> ids(N);
    std::vector<real> pos_x(N), pos_y(N), pos_z(N);
    std::vector<real> vel_x(N), vel_y(N), vel_z(N);
    std::vector<real> acc_x(N), acc_y(N), acc_z(N);
    std::vector<real> mass(N), dens(N), pres(N), ene(N), sml(N), sound(N);
    std::vector<real> alpha(N), balsara(N), gradh(N), phi(N);
    std::vector<int> neighbor(N);
    
    for (hsize_t i = 0; i < N; ++i) {
        const SPHParticle* p = particles[i];
        ids[i] = p->id;
        pos_x[i] = p->pos[0];
        pos_y[i] = p->pos[1];
        pos_z[i] = p->pos[2];
        vel_x[i] = p->vel[0];
        vel_y[i] = p->vel[1];
        vel_z[i] = p->vel[2];
        acc_x[i] = p->acc[0];
        acc_y[i] = p->acc[1];
        acc_z[i] = p->acc[2];
        mass[i] = p->mass;
        dens[i] = p->dens;
        pres[i] = p->pres;
        ene[i] = p->ene;
        sml[i] = p->sml;
        sound[i] = p->sound;
        alpha[i] = p->alpha;
        balsara[i] = p->balsara;
        gradh[i] = p->gradh;
        phi[i] = p->phi;
        neighbor[i] = p->neighbor;
    }
    
    // Write datasets
    write_1d_int_dataset("id", ids);
    write_1d_dataset("pos_x", pos_x);
    write_1d_dataset("pos_y", pos_y);
    write_1d_dataset("pos_z", pos_z);
    write_1d_dataset("vel_x", vel_x);
    write_1d_dataset("vel_y", vel_y);
    write_1d_dataset("vel_z", vel_z);
    write_1d_dataset("acc_x", acc_x);
    write_1d_dataset("acc_y", acc_y);
    write_1d_dataset("acc_z", acc_z);
    write_1d_dataset("mass", mass);
    write_1d_dataset("dens", dens);
    write_1d_dataset("pres", pres);
    write_1d_dataset("ene", ene);
    write_1d_dataset("sml", sml);
    write_1d_dataset("sound", sound);
    write_1d_dataset("alpha", alpha);
    write_1d_dataset("balsara", balsara);
    write_1d_dataset("gradh", gradh);
    write_1d_dataset("phi", phi);
    write_1d_int_dataset("neighbor", neighbor);
    
    if (plist_id != H5P_DEFAULT) {
        H5Pclose(plist_id);
    }
    H5Gclose(group_id);
    return true;
}

bool HDF5Writer::write_energy_group(const OutputMetadata& metadata) {
    // Create energy group
    hid_t group_id = H5Gcreate(m_file_id, "/energy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0) return false;
    
    // Write scalar datasets
    hsize_t dims[1] = {1};
    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    
    hid_t dataset_id;
    
    dataset_id = H5Dcreate(group_id, "kinetic", H5T_NATIVE_DOUBLE, dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &metadata.kinetic_energy);
    H5Dclose(dataset_id);
    
    dataset_id = H5Dcreate(group_id, "thermal", H5T_NATIVE_DOUBLE, dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &metadata.thermal_energy);
    H5Dclose(dataset_id);
    
    dataset_id = H5Dcreate(group_id, "potential", H5T_NATIVE_DOUBLE, dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &metadata.potential_energy);
    H5Dclose(dataset_id);
    
    dataset_id = H5Dcreate(group_id, "total", H5T_NATIVE_DOUBLE, dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &metadata.total_energy);
    H5Dclose(dataset_id);
    
    H5Sclose(dataspace_id);
    H5Gclose(group_id);
    return true;
}

bool HDF5Writer::write_string_attribute(hid_t loc_id, const char* name, const std::string& value) {
    hid_t dataspace_id = H5Screate(H5S_SCALAR);
    hid_t string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type, value.size() + 1);
    H5Tset_strpad(string_type, H5T_STR_NULLTERM);
    
    hid_t attr_id = H5Acreate(loc_id, name, string_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id < 0) {
        H5Tclose(string_type);
        H5Sclose(dataspace_id);
        return false;
    }
    
    H5Awrite(attr_id, string_type, value.c_str());
    
    H5Aclose(attr_id);
    H5Tclose(string_type);
    H5Sclose(dataspace_id);
    return true;
}

bool HDF5Writer::read_string_attribute(hid_t loc_id, const char* name, std::string& value) {
    if (H5Aexists(loc_id, name) <= 0) {
        return false;
    }
    
    hid_t attr_id = H5Aopen(loc_id, name, H5P_DEFAULT);
    if (attr_id < 0) return false;
    
    hid_t string_type = H5Aget_type(attr_id);
    size_t str_size = H5Tget_size(string_type);
    
    std::vector<char> buffer(str_size + 1);
    H5Aread(attr_id, string_type, buffer.data());
    value = std::string(buffer.data());
    
    H5Tclose(string_type);
    H5Aclose(attr_id);
    return true;
}

bool HDF5Writer::read_metadata(const std::string& filepath, OutputMetadata& metadata) {
    // Open file for reading
    hid_t file_id = H5Fopen(filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) return false;
    
    // Open metadata group
    hid_t group_id = H5Gopen(file_id, "/metadata", H5P_DEFAULT);
    if (group_id < 0) {
        H5Fclose(file_id);
        return false;
    }
    
    // Read JSON metadata if available (easiest way)
    std::string metadata_json_str;
    if (read_string_attribute(group_id, "metadata_json", metadata_json_str)) {
        nlohmann::json j = nlohmann::json::parse(metadata_json_str);
        metadata = OutputMetadata::from_json(j);
    }
    
    H5Gclose(group_id);
    H5Fclose(file_id);
    return true;
}

bool HDF5Writer::read_particles(const std::string& filepath, std::vector<SPHParticle*>& particles) {
    // Open file for reading
    hid_t file_id = H5Fopen(filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) return false;
    
    // Open particles group
    hid_t group_id = H5Gopen(file_id, "/particles", H5P_DEFAULT);
    if (group_id < 0) {
        H5Fclose(file_id);
        return false;
    }
    
    // Get particle count from id dataset
    hid_t dataset_id = H5Dopen(group_id, "id", H5P_DEFAULT);
    if (dataset_id < 0) {
        H5Gclose(group_id);
        H5Fclose(file_id);
        return false;
    }
    
    hid_t dataspace_id = H5Dget_space(dataset_id);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
    const hsize_t N = dims[0];
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    
    // Helper lambda to read 1D dataset
    auto read_1d_dataset = [&](const char* name, std::vector<real>& data) -> bool {
        data.resize(N);
        hid_t dset_id = H5Dopen(group_id, name, H5P_DEFAULT);
        if (dset_id < 0) return false;
        H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(dset_id);
        return true;
    };
    
    // Helper lambda to read 1D integer dataset
    auto read_1d_int_dataset = [&](const char* name, std::vector<int>& data) -> bool {
        data.resize(N);
        hid_t dset_id = H5Dopen(group_id, name, H5P_DEFAULT);
        if (dset_id < 0) return false;
        H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(dset_id);
        return true;
    };
    
    // Read all datasets
    std::vector<int> ids(N);
    std::vector<real> pos_x(N), pos_y(N), pos_z(N);
    std::vector<real> vel_x(N), vel_y(N), vel_z(N);
    std::vector<real> acc_x(N), acc_y(N), acc_z(N);
    std::vector<real> mass(N), dens(N), pres(N), ene(N), sml(N), sound(N);
    std::vector<real> alpha(N), balsara(N), gradh(N), phi(N);
    std::vector<int> neighbor(N);
    
    read_1d_int_dataset("id", ids);
    read_1d_dataset("pos_x", pos_x);
    read_1d_dataset("pos_y", pos_y);
    read_1d_dataset("pos_z", pos_z);
    read_1d_dataset("vel_x", vel_x);
    read_1d_dataset("vel_y", vel_y);
    read_1d_dataset("vel_z", vel_z);
    read_1d_dataset("acc_x", acc_x);
    read_1d_dataset("acc_y", acc_y);
    read_1d_dataset("acc_z", acc_z);
    read_1d_dataset("mass", mass);
    read_1d_dataset("dens", dens);
    read_1d_dataset("pres", pres);
    read_1d_dataset("ene", ene);
    read_1d_dataset("sml", sml);
    read_1d_dataset("sound", sound);
    read_1d_dataset("alpha", alpha);
    read_1d_dataset("balsara", balsara);
    read_1d_dataset("gradh", gradh);
    read_1d_dataset("phi", phi);
    read_1d_int_dataset("neighbor", neighbor);
    
    // Clear existing particles and create new ones
    for (auto* p : particles) {
        delete p;
    }
    particles.clear();
    particles.reserve(N);
    
    for (hsize_t i = 0; i < N; ++i) {
        SPHParticle* p = new SPHParticle();
        p->id = ids[i];
        p->pos[0] = pos_x[i];
        p->pos[1] = pos_y[i];
        p->pos[2] = pos_z[i];
        p->vel[0] = vel_x[i];
        p->vel[1] = vel_y[i];
        p->vel[2] = vel_z[i];
        p->acc[0] = acc_x[i];
        p->acc[1] = acc_y[i];
        p->acc[2] = acc_z[i];
        p->mass = mass[i];
        p->dens = dens[i];
        p->pres = pres[i];
        p->ene = ene[i];
        p->sml = sml[i];
        p->sound = sound[i];
        p->alpha = alpha[i];
        p->balsara = balsara[i];
        p->gradh = gradh[i];
        p->phi = phi[i];
        p->neighbor = neighbor[i];
        p->next = nullptr;
        
        particles.push_back(p);
    }
    
    H5Gclose(group_id);
    H5Fclose(file_id);
    return true;
}

} // namespace sph
