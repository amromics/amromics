/* eslint-disable */
import axios from 'axios'


export default {
  fetchResult (sampleId) {
        console.log(process.env.BASE_URL)
        return axios.get('/static/data/samples/' + sampleId + '.json')
  },
  fetchSetResult () {
        console.log(process.env.BASE_URL)
        return axios.get('/static/data/set.json')
  }
}
