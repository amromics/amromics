/* eslint-disable */
import axios from 'axios'


export default {
    fetchResult(collectionID, sampleId) {
        console.log(process.env.BASE_URL)
        return axios.get('../../static/data/' + collectionID + '/samples/' + sampleId + '.json')
    },
    fetchSetResult(collectionID) {
        return axios.get('../static/data/' + collectionID + '/set.json')
    },
    fetchListCollectionResult() {
        return axios.get('../static/data/collections.json')
    },
    fetchPangenomeCluster(collectionID) {
        return axios.get('../static/data/' + collectionID + '/set/pangenome_cluster.json')
    },
    fetchPangenomeSummary(collectionID) {
        return axios.get('../static/data/' + collectionID + '/set/pangenome_summary.json')
    },
    fetchHeatmap(collectionID) {
        return axios.get('../static/data/' + collectionID + '/set/amrheatmap.json')
    },
    async fetchAlignment(collectionID, gene) {

        // axios.get('/static/data/' + collectionID + '/set/alignments/' + gene + '.json.gz', { responseType: 'arraybuffer' }).then(function (response) {
        //     //console.log(pako.ungzip(response.data, { to: 'string' }))
        //     const zlib = require("zlib");
        //     zlib.gunzip(response.data, function(_err, output) {
        //         console.log(_err.toString());
        //         console.log(output.toString());
        //     });
        // });;
       // const pako = require('pako');
        return axios({
            method:'get',
            url:'../static/data/' + collectionID + '/set/alignments/' + gene + '.json.gz',
            responseType:'arraybuffer'
            })
           
    }
}
